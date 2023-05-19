from __future__ import annotations

import datetime
import argparse
import logging
import uuid
from requests import post
from requests.auth import HTTPBasicAuth

import sbol3
import tyto
from sbol_utilities.workarounds import type_to_standard_extension

COMPLEXITY_SCORE_NAMESPACE = 'http://igem.org/IDT_complexity_score'
REPORT_ACTIVITY_TYPE = 'https://github.com/SynBioDex/SBOL-utilities/compute-sequence-complexity'


class IDTAccountAccessor:
    """Class that wraps access to the IDT API"""

    _TOKEN_URL = 'https://www.idtdna.com/Identityserver/connect/token'
    """API URL for obtaining session tokens"""
    _SCORE_URL = 'https://www.idtdna.com/api/complexities/screengBlockSequences'
    """APR URL for obtaining sequence scores"""
    _BLOCK_SIZE = 1  # TODO: determine if it is possible to run multiple sequences in a single query
    SCORE_TIMEOUT = 120
    """Number of seconds to wait for score query requests to complete"""

    def __init__(self, username: str, password: str, client_id: str, client_secret: str):
        """Initialize with required access information for IDT API (see: https://www.idtdna.com/pages/tools/apidoc)
        Automatically logs in and obtains a session token

        :param username: Username of your IDT account
        :param password: Password of your IDT account
        :param client_id: ClientID key of your IDT account
        :param client_secret: ClientSecret key of your IDT account
        """
        self.username = username
        self.password = password
        self.client_id = client_id
        self.client_secret = client_secret
        self.token = self._get_idt_access_token()

    @staticmethod
    def from_json(json_object) -> IDTAccountAccessor:
        """Initialize IDT account accessor from a JSON object with field values

        :param json_object: object with account information
        :return: Account accessor object
        """
        return IDTAccountAccessor(username=json_object['username'], password=json_object['password'],
                                  client_id=json_object['ClientID'], client_secret=json_object['ClientSecret'])

    def _get_idt_access_token(self) -> str:
        """Get access token for IDT API (see: https://www.idtdna.com/pages/tools/apidoc)

        :return: access token string
        """
        logging.info('Connecting to IDT API')
        data = {'grant_type': 'password', 'username': self.username, 'password': self.password, 'scope': 'test'}
        auth = HTTPBasicAuth(self.client_id, self.client_secret)
        result = post(IDTAccountAccessor._TOKEN_URL, data, auth=auth, timeout=IDTAccountAccessor.SCORE_TIMEOUT)

        if 'access_token' in result.json():
            return result.json()['access_token']
        else:
            raise ValueError('Access token for IDT API could not be generated. Check your credentials.')

    def get_sequence_scores(self, sequences: list[sbol3.Sequence]) -> list:
        """Retrieve synthesis complexity scores of sequences from the IDT API
        This system uses the gBlock API, which is intended for sequences from 125 to 3000 bp in length

        :param sequences: sequences for which we want to calculate the complexity score
        :return: dictionary mapping sequences to complexity Scores
        :return: List of lists of dictionaries with information about sequence synthesizability features
        """
        # Set up list of query dictionaries
        seq_dict = [{'Name': str(seq.display_name), 'Sequence': str(seq.elements)} for seq in sequences]
        # Break into query blocks
        partitions_sequences = [seq_dict[x:x + 1] for x in range(0, len(seq_dict), IDTAccountAccessor._BLOCK_SIZE)]
        # Send each query to IDT and collect results
        results = []
        for idx, partition in enumerate(partitions_sequences):
            logging.debug('Sequence score request %i of %i', idx+1, len(partitions_sequences))
            resp = post(IDTAccountAccessor._SCORE_URL, json=partition, timeout=IDTAccountAccessor.SCORE_TIMEOUT,
                        headers={'Authorization': 'Bearer {}'.format(self.token),
                                 'Content-Type': 'application/json; charset=utf-8'})
            response_list = resp.json()
            if len(response_list) != len(partition):
                raise ValueError(f'Unexpected complexity score: expected {len(partition)} scores, '
                                 f'but got {len(response_list)}')
            results.append(resp.json())
        logging.info('Requests to IDT API finished.')
        return results

    def get_sequence_complexity(self, sequences: list[sbol3.Sequence]) -> dict[sbol3.Sequence, float]:
        """ Extract complexity scores from IDT API for a list of SBOL Sequence objects
        This works by computing full sequence evaluations, then compressing down to a single score for each sequence.

        :param sequences: list of SBOL Sequences to evaluate
        :return: dictionary mapping sequences to complexity Scores
        """
        # Retrieve full evaluations for sequences
        scores = self.get_sequence_scores(sequences)
        # Compute total score for each sequence as the sum all complexity scores for the sequence
        score_list = []
        for score_set in scores:
            for sequence_scores in score_set:
                complexity_score = sum(score.get('Score') for score in sequence_scores)
                score_list.append(complexity_score)
        # Associate each sequence to its score
        return dict(zip(sequences, score_list))

# TODO: separate scoring sequences from scoring documents
# TODO: function for retrieving sequence scores


def IDT_calculate_complexity_score(username: str, password: str, client_id: str, client_secret: str,
                                   doc: sbol3.Document) -> dict[sbol3.Sequence, float]:
    """ Add computed complexity scores of sequences to SBOL document with a timestamp

    :param username: Username of your IDT account
    :param password: Password of your IDT account
    :param client_id: ClientID key of your IDT account
    :param client_secret: ClientSecret key of your IDT account
    :param doc: SBOL document with sequences of interest in it
    :return: Dictionary mapping Sequences to complexity scores
    """
    # Extract sequence identities and elements from SBOL document
    logging.info(f'Importing distribution sequences')
    sequences = [top_level for top_level in doc if isinstance(top_level, sbol3.Sequence)]
    # TODO: only compute for sequences that don't already have scores
    # TODO: separate doc scrape from sequence score computation
    # TODO: wrap IDT account into info an object
    # Query for the scores of the sequences
    idt_accessor = IDTAccountAccessor(username, password, client_id, client_secret)
    score_dictionary = idt_accessor.get_sequence_complexity(sequences)

    # Create report generation activity
    timestamp = datetime.datetime.utcnow().isoformat(timespec='seconds') + 'Z'
    report_id = f'{COMPLEXITY_SCORE_NAMESPACE}/Complexity_Report_{timestamp.replace(":", "").replace("-", "")}_' \
                f'{str(uuid.uuid4())[0:8]}'
    report_generation = sbol3.Activity(report_id, end_time=timestamp, types=[REPORT_ACTIVITY_TYPE])
    doc.add(report_generation)
    # Mark the sequences with their scores, where each score is a dimensionless measure
    for sequence, score in score_dictionary.items():
        measure = sbol3.Measure(score, unit=tyto.OM.number_unit, types=[tyto.EDAM.sequence_complexity_report])
        measure.generated_by.append(report_generation)
        sequence.measures.append(measure)

    return score_dictionary


def main():
    """
    Main wrapper: read from input file, invoke IDT_calculate_complexity_score, then write to output file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('username', help="Username of your IDT account")
    parser.add_argument('password', help="Password of your IDT account")
    parser.add_argument('ClientID', help="ClientID of your IDT account")
    parser.add_argument('ClientSecret', help="ClientSecret of your IDT account")
    parser.add_argument('input_file', help="Absolute path to sbol file with sequences")
    parser.add_argument('output_name', help="Name of SBOL file to be written")
    parser.add_argument('-t', '--file-type', dest='file_type', default=sbol3.SORTED_NTRIPLES,
                        help="Name of SBOL file to output to (excluding type)")
    parser.add_argument('--verbose', '-v', dest='verbose', action='count', default=0)
    args_dict = vars(parser.parse_args())

    # Extract arguments:
    verbosity = args_dict['verbose']
    logging.getLogger().setLevel(level=(logging.WARN if verbosity == 0 else
                                        logging.INFO if verbosity == 1 else logging.DEBUG))
    input_file = args_dict['input_file']
    output_name = args_dict['output_name']

    extension = type_to_standard_extension[args_dict['file_type']]
    outfile_name = output_name if output_name.endswith(extension) else output_name + extension

    # Read file, convert, and write resulting document
    logging.info('Reading SBOL file ' + input_file)
    doc = sbol3.Document()
    doc.read(input_file)
    results = IDT_calculate_complexity_score(args_dict['username'], args_dict['password'], args_dict['ClientID'], args_dict['ClientSecret'], doc)
    doc.write(outfile_name, args_dict['file_type'])
    logging.info('SBOL file written to %s with %i new scores calculated', outfile_name, len(results))


if __name__ == '__main__':
    main()
