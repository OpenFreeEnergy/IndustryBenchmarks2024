# Preperation of systems

Done by J. Bluck

<ol type="1">
<li>Read in original inputs from the repo
<li>Duplicate the input, to ensure there is a reference to comapre to
<li>If multiple ligands, extract preferred ligand input (after discussion with the OpenFE team)
<li>Extract all cofactors into a single maestro entry
<li>Check the sequence to see if it contains any non-natural amino acids, mutate back if I do
<li>Check all caps to see if any need to be charged or not. Leave residue if charged
<li>Open the Protein preparation wizard, if caps are needed
<li>Run only step 3 - "Preprocess". The only three options activated are:
   <ul type="disc">
    <li>Cap termini
    <li>Convert selomethionines to methionines
    <li>Include peptides when capping termini
    </ul>
<li>Change the cap names if not ACE, NME:
   <ul type="disc">
    <li>Highlight all waters
    <li>Builder > Other Edits > Change Atom Properties > Property: Residue/Chain Name
    <li>Residue name "NME"
    </ul>
<li>Export outputs, named "ligand.sdf", "protein.pdb" and optionally "cofactors.sdf"
<li>Run input through validation script
<li>Scrub files of sensitive information
</ol>

NOTE: There is a crystallographic GOL residue that has been isolated as a cofactor. But this should probably be excldued from future calculations.
