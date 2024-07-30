# Preperation of systems

Done by J. Bluck

<ol type="1">
<li>Read in original inputs from the repo
<li>Duplicate the input, to ensure there is a reference to comapre to
<li>Check the sequence to see if it contains any non-natural amino acids, mutate back if I do
<li>Open the Protein preparation wizard (Schrodinger 2024-2), if caps are needed
<li>Run only step 3 - "Preprocess". The only three options activated are:
   <ul type="disc">
    <li>Cap termini
    <li>Convert selomethionines to methionines
    <li>Include peptides when capping termini
    </ul>
<li>Convert Ser terminal back to charged
<li>Change the cap names if not ACE, NME:
   <ul type="disc">
    <li>Highlight all waters
    <li>Builder > Other Edits > Change Atom Properties > Property: Residue/Chain Name
    <li>Residue name "NME"
    </ul>
<li>Export outputs, named "ligand.sdf", "protein.pdb"
<li>Run input through validation script
<li>Scrub files of sensitive information
</ol>
