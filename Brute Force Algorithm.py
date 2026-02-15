def get_amino_acid_masses():
    """Returns dictionary of amino acid masses"""
    return {
        'G':57, 'A':71, 'S':87, 'P':97, 'V':99,
        'T':101, 'C':103, 'I':113, 'L':113, 'N':114,
        'D':115, 'K':128, 'Q':128, 'E':129, 'M':131,
        'H':137, 'F':147, 'R':156, 'Y':163, 'W':186
    }

def get_mass_to_amino():
    """Returns dictionary mapping masses to amino acids"""
    aa_masses=get_amino_acid_masses()
    mass_to_aa={}
    for aa,mass in aa_masses.items():
        if mass not in mass_to_aa:
            mass_to_aa[mass]=[]
        mass_to_aa[mass].append(aa)
    return mass_to_aa

def linear_spectrum(peptide,aa_masses):
    """Generate theoretical linear spectrum of a peptide"""
    prefix_mass=[0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1]+aa_masses[aa])
    linear_spec=[0]
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            linear_spec.append(prefix_mass[j]-prefix_mass[i])
    return sorted(linear_spec)

def cyclic_spectrum(peptide,aa_masses):
    """Generate theoretical cyclic spectrum of a peptide"""
    prefix_mass=[0]
    for aa in peptide:
        prefix_mass.append(prefix_mass[-1]+aa_masses[aa])
    peptide_mass=prefix_mass[-1]
    cyclic_spec=[0]
    for i in range(len(peptide)):
        for j in range(i+1,len(peptide)+1):
            cyclic_spec.append(prefix_mass[j]-prefix_mass[i])
            if i>0 and j<len(peptide):
                cyclic_spec.append(peptide_mass-(prefix_mass[j]-prefix_mass[i]))
    return sorted(cyclic_spec)

def is_consistent(peptide,spectrum,aa_masses):
    """Check if linear spectrum of peptide is consistent with given spectrum"""
    lin_spec=linear_spectrum(peptide,aa_masses)
    #Count occurrences in both spectra
    spec_count={}
    for mass in spectrum:
        spec_count[mass]=spec_count.get(mass,0)+1
    for mass in lin_spec:
        if mass not in spec_count or spec_count[mass]==0:
            return False
        spec_count[mass]-=1
    return True

def peptide_mass(peptide,aa_masses):
    """Calculate total mass of peptide"""
    return sum(aa_masses[aa] for aa in peptide)

def brute_force_cyclopeptide_sequencing(spectrum):
    """Brute Force Cyclopeptide Sequencing Algorithm"""
    aa_masses=get_amino_acid_masses()
    mass_to_aa=get_mass_to_amino()
    #Get unique amino acid masses that appear in the spectrum
    possible_masses=set()
    for mass in spectrum:
        if mass in mass_to_aa:
            possible_masses.add(mass)
    print("="*80)
    print("BRUTE FORCE CYCLOPEPTIDE SEQUENCING")
    print("="*80)
    print(f"\nGiven Spectrum: {spectrum}")
    print(f"Parent Mass: {max(spectrum)}")
    print(f"\nPossible amino acid masses in spectrum: {sorted(possible_masses)}")
    print("Corresponding amino acids:")
    for mass in sorted(possible_masses):
        print(f"  Mass {mass}: {', '.join(mass_to_aa[mass])}")
    parent_mass=max(spectrum)
    candidate_peptides=[''] #Start with empty peptide
    final_peptides=[]
    iteration=0
    while candidate_peptides:
        iteration+=1
        print(f"\n{'='*80}")
        print(f"ITERATION {iteration}")
        print(f"{'='*80}")
        print(f"Number of candidate peptides to extend: {len(candidate_peptides)}")
        print(f"Current candidates: {candidate_peptides if candidate_peptides!=[''] else ['(empty)']}")
        new_candidates=[]
        #Step 1: Expand all candidate peptides by one amino acid
        print(f"\n--- STEP 1: EXTENDING PEPTIDES ---")
        print(f"Extending each candidate with amino acids: P(97), V(99), T(101), C(103)")
        print()
        for peptide in candidate_peptides:
            display_peptide=peptide if peptide else "(empty)"
            print(f"Extending '{display_peptide}':")
            for mass in sorted(possible_masses):
                for aa in mass_to_aa[mass]:
                    new_peptide=peptide+aa
                    new_candidates.append(new_peptide)
                    new_mass=peptide_mass(new_peptide,aa_masses)
                    print(f"  + {aa}({mass}) → '{new_peptide}' [mass={new_mass}]")
            print()
        print(f"Total extended peptides generated: {len(new_candidates)}")
        print(f"All extended peptides: {new_candidates}")
        #Step 2: Form theoretical spectra and check consistency
        print(f"\n--- STEP 2: CHECKING CONSISTENCY WITH SPECTRUM ---")
        candidate_peptides=[]
        for peptide in new_candidates:
            current_mass=peptide_mass(peptide,aa_masses)
            if current_mass==parent_mass:
                #Check if cyclic spectrum matches
                cyc_spec=cyclic_spectrum(peptide,aa_masses)
                print(f"\n'{peptide}' [mass={current_mass}] - COMPLETE PEPTIDE")
                print(f"  Cyclic spectrum: {cyc_spec}")
                if cyc_spec==spectrum:
                    final_peptides.append(peptide)
                    print(f"  ✓✓✓ MATCH! This peptide matches the given spectrum! ✓✓✓")
                else:
                    print(f"  ✗ No match - spectrum doesn't match")
            elif current_mass<parent_mass:
                #Check if linear spectrum is consistent
                lin_spec=linear_spectrum(peptide,aa_masses)
                consistent=is_consistent(peptide,spectrum,aa_masses)
                print(f"\n'{peptide}' [mass={current_mass}]")
                print(f"  Linear spectrum: {lin_spec}")
                if consistent:
                    candidate_peptides.append(peptide)
                    print(f"  ✓ Consistent - keeping for next iteration")
                else:
                    print(f"  ✗ Not consistent - discarding")
            else:
                #Mass exceeds parent mass
                print(f"\n'{peptide}' [mass={current_mass}]")
                print(f"  ✗ Mass exceeds parent mass - discarding")
        print(f"\n--- END OF ITERATION {iteration} ---")
        print(f"Remaining candidates for next iteration: {len(candidate_peptides)}")
        if len(candidate_peptides)>0:
            print(f"Candidates: {candidate_peptides}")
        print(f"Complete matches found so far: {len(final_peptides)}")
        if final_peptides:
            print(f"Matched peptides: {final_peptides}")
    print(f"\n{'='*80}")
    print("FINAL RESULTS")
    print(f"{'='*80}")
    print(f"\nTotal peptides found: {len(final_peptides)}")
    if final_peptides:
        print("\n✓✓✓ FINAL PEPTIDE SEQUENCES ✓✓✓")
        for i,peptide in enumerate(final_peptides,1):
            mass=peptide_mass(peptide,aa_masses)
            cyc_spec=cyclic_spectrum(peptide,aa_masses)
            print(f"\n{i}. Peptide: {peptide}")
            print(f"   Mass: {mass}")
            print(f"   Cyclic Spectrum: {cyc_spec}")
    else:
        print("\n✗ No peptides found matching the spectrum!")
    return final_peptides

#Given spectrum from the slide
spectrum=[0,97,97,99,101,103,196,198,198,200,202,295,297,299,299,301,394,396,398,400,400,497]

print("\nRunning Brute Force Cyclopeptide Sequencing Algorithm")
print("="*80)

result=brute_force_cyclopeptide_sequencing(spectrum)

print("\n"+"="*80)
print("ALGORITHM COMPLETE")
print("="*80)
if result:
    print(f"Found {len(result)} matching peptide(s): {result}")
