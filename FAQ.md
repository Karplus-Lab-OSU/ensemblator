FAQ:

	Q:	Do they have to be NMR ensembles?
	A:	No. You can make an ensemble out of any type of structure you want!
	
	Q:	My pdb files keep throwing errors after refinement with phenix... what's the deal?
	A:	Phenix loves to introduce lowercase letters into pdb files, if you have any heavy atoms
		in your structure (ie. Ca, Cl, Br, etc.). These cause biopython to freak out, and will
		break the ensemblator. Try making everything in the file uppercase, or try eliminating
		all the non-ATOM lines from the pdb file.
