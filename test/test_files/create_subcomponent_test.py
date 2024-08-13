from sbol3 import *
from sbol_utilities.calculate_sequences import compute_sequence
from sbol_utilities.component import *
from sbol_utilities.helper_functions import url_to_identity
import tyto

set_namespace('https://synbiohub.org/public/igem/')
doc = Document()

i13504 = Component('i13504', SBO_DNA)
i13504.name = 'iGEM 2016 interlab reporter'
i13504.description = 'GFP expression cassette used for 2016 iGEM interlab study'
i13504.roles.append(tyto.SO.engineered_region)

doc.add(i13504)

b0034, b0034_seq = doc.add(rbs('B0034', sequence='aaagaggagaaa', name='RBS (Elowitz 1999)'))


e0040_sequence = 'atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaataataa'
e0040, _ = doc.add(cds('E0040', sequence=e0040_sequence, name='GFP'))

b0015_sequence = 'ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctctactagagtcacactggctcaccttcgggtgggcctttctgcgtttata'
b0015, _ = doc.add(terminator('B0015', sequence=b0015_sequence, name='double terminator'))


order(b0034, e0040, i13504)
order(e0040, b0015, i13504)

i13504_seq = compute_sequence(i13504)

b0015_subcomponent = next(f for f in i13504.features if f.instance_of == b0015.identity)
b0015_range = b0015_subcomponent.locations[0]
print(f'Range of {b0015.display_name}: ({b0015_range.start}, {b0015_range.end})')

i13504_system = functional_component('i13504_system')
doc.add(i13504_system)

gfp = add_feature(i13504_system, ed_protein('https://www.fpbase.org/protein/gfpmut3/', name='GFP'))

i13504_subcomponent = add_feature(i13504_system, i13504)

e0040_subcomponent = next(f for f in i13504.features if f.instance_of == e0040.identity)
e0040_reference = ComponentReference(i13504_subcomponent, e0040_subcomponent)
i13504_system.features.append(e0040_reference)


add_interaction(tyto.SBO.genetic_production,
                participants={gfp: tyto.SBO.product, e0040_reference: tyto.SBO.template})

J23101_sequence = 'tttacagctagctcagtcctaggtattatgctagc'
J23101, _ = doc.add(promoter('J23101', sequence=J23101_sequence))
J23106_sequence = 'tttacggctagctcagtcctaggtatagtgctagc'
J23106, _ = doc.add(promoter('J23106', sequence=J23106_sequence))

device1 = doc.add(functional_component('interlab16device1'))
device1_i13504_system = add_feature(device1, SubComponent(i13504_system))
order(J23101, ComponentReference(device1_i13504_system, i13504_subcomponent), device1)
device2 = doc.add(functional_component('interlab16device2'))
device2_i13504_system = add_feature(device2, SubComponent(i13504_system))
order(J23106, ComponentReference(device2_i13504_system, i13504_subcomponent), device2)
print(f'Device 1 second subcomponent points to {device1.constraints[0].object.lookup().refers_to.lookup().instance_of}')

interlab16 = doc.add(Collection('interlab16',members=[device1, device2]))
print(f'Members are {", ".join(m.lookup().display_id for m in interlab16.members)}')

ecoli = doc.add(strain('Ecoli_DH5_alpha'))
pSB1C3 = doc.add(Component('pSB1C3', SBO_DNA, roles=[tyto.SO.plasmid_vector]))

device1_ecoli = doc.add(strain('device1_ecoli'))

plasmid = LocalSubComponent(SBO_DNA, roles=[tyto.SO.plasmid_vector], name="Interlab Device 1 in pSB1C3")
device1_ecoli.features.append(plasmid)
device1_subcomponent = contains(plasmid, device1)
contains(plasmid, pSB1C3)
order(device1, pSB1C3, device1_ecoli)

contains(ecoli, plasmid, device1_ecoli)


report = doc.validate()
if report:
    print('Document is not valid')
    print(f'Document has {len(report.errors)} errors')
    print(f'Document has {len(report.warnings)} warnings')
else:
    print('Document is valid')

doc.write('subcomponent_test.nt', file_format=SORTED_NTRIPLES)
