
import os, sys, re, csv, glob
import Bio.PDB
import numpy as np


usage = '\nUsage: python Filter_BindCraft_designs.py bindcraft_folder_path min_iptm\nRemoves any with iptm < min_iptm, keeps only one seq per traj, checks secstruct, and checks termini\n'


#NOTE: DESIGNS SHOULD ALREADY BE RANKED BY IPTM IN final_design_stats.csv TO ENSURE BEST ONE FOR EACH TRAJ GETS KEPT


bindcraft_folder_path = sys.argv[1]
min_iptm = float(sys.argv[2])

with open(os.path.join(bindcraft_folder_path, 'final_design_stats.csv'), 'r') as csvfile:
	csv_lines = [line for line in csvfile]


csv_lines = csv_lines[1:]

csv_lines = csv.reader(csv_lines, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
csv_lines = [l for l in csv_lines]


parser = Bio.PDB.PDBParser(QUIET = True)


#in chain B
def count_helices(path_to_pdb):
	structure = parser.get_structure('scratch', path_to_pdb)
	model = structure[0]
	dssp = Bio.PDB.DSSP(model, path_to_pdb)
	secstruct = ''
	for i in range(len(model['A']), len(model['A']) + len(model['B'])): #from https://github.com/lamm-mit/ProteinDiffusionGenerator/blob/main/ProteinDiffusionGenerator_Model_B.ipynb
		a_key = list(dssp.keys())[i]
		secstruct += dssp[a_key][2]
	helices = re.findall('H' + '+', secstruct)
	sheets = re.findall('E' + '+', secstruct)
	return len(helices), len(sheets)


class Binder:
	def __init__(self, csv_line):
		self.csv_line = csv_line
		self.iptm = float(csv_line[23])
		self.ipae = float(csv_line[35])
		self.pdbpath = glob.glob(os.path.join(bindcraft_folder_path, 'Accepted/Ranked', '*' + csv_line[1] + '*.pdb'))[0]
		self.traj = csv_line[1].split('_mpnn')[0]
		self.name = csv_line[1]

	def secstruct_check(self):
		numhlx, numsht = count_helices(self.pdbpath)
		if (numsht >= 3) or (numhlx >= 3):
			return True
		else:
			return False


binders = [Binder(csvl) for csvl in csv_lines]

binders_iptm = [binder for binder in binders if binder.iptm >= min_iptm]


unique_trajectories = []
binders_iptm_oneseqpertraj = []

for binder in binders_iptm:
	if binder.traj not in unique_trajectories:
		binders_iptm_oneseqpertraj.append(binder)
		unique_trajectories.append(binder.traj)


binders_iptm_oneseqpertraj_secstruct = [binder for binder in binders_iptm_oneseqpertraj if binder.secstruct_check()]


#Next do termini checking...

#Credit for contact map: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two) :
	"""Returns a matrix of C-alpha distances between two chains"""
	answer = np.zeros((len(chain_one), len(chain_two)), float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			answer[row, col] = calc_residue_dist(residue_one, residue_two)
	return answer


def get_A_B_contact_map(path_to_pdb, cutoff_distance):
	structure = parser.get_structure('scratch', path_to_pdb)
	model = structure[0]
	dist_matrix = calc_dist_matrix(model['A'], model['B'])
	contact_map = dist_matrix < cutoff_distance
	return contact_map


#n from 0, ctmapT = transpose of ctmap
def is_res_n_of_B_in_contact_with_A(ctmapT, n):
	if True in ctmapT[n]:
		return True
	else:
		return False


contact_cutoff_distance = 10.0


def Nterm_Cterm_free(path_to_pdb):
	ctm = get_A_B_contact_map(path_to_pdb, contact_cutoff_distance)
	ctmT = np.transpose(ctm) #Now rows = res in chain B
	nterm_free = True
	cterm_free = True
	if True in ctmT[0]:
		nterm_free = False
	if True in ctmT[-1]:
		cterm_free = False
	return nterm_free, cterm_free


binders_iptm_oneseqpertraj_secstruct_ntermfree = []
termini_csv_lines = ['Design,avg iptm,avg ipae,Nterm_free,Cterm_free\n']

filpath = os.path.join(bindcraft_folder_path, 'filtered')
os.system(f'mkdir {filpath}')

for binder in binders_iptm_oneseqpertraj_secstruct:
	ntf, ctf = Nterm_Cterm_free(binder.pdbpath)
	if ntf:
		binders_iptm_oneseqpertraj_secstruct_ntermfree.append(binder)
		termini_csv_lines.append(f'{binder.name},{binder.iptm},{binder.ipae},{ntf},{ctf}\n')
		os.system(f'cp {binder.pdbpath} {filpath}')


txt_header = [f'Filtered {bindcraft_folder_path} with min_iptm {min_iptm}\n', f'{bindcraft_folder_path} had {len(binders)} accepted binders.\n', f'{len(binders_iptm)} had iptm >= {min_iptm}.\n', f'Keeping one seq per traj leaves {len(binders_iptm_oneseqpertraj)}.\n', f'{len(binders_iptm_oneseqpertraj_secstruct)} of those met secstruct requirements.\n', f'{len(binders_iptm_oneseqpertraj_secstruct_ntermfree)} of those had the N-terminus free.\n']

termini_csv_lines = txt_header + ['\n'] + termini_csv_lines

with open(os.path.join(filpath, 'filtered.txt'), 'w') as txt:
	for line in termini_csv_lines:
		txt.write(line)

print()
for line in txt_header:
	print(line.strip())
print()


