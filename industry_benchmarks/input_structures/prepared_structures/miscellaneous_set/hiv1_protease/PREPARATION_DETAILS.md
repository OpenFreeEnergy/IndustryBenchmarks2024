# Preperation of systems

Done by J. Bluck

<ol type="1">
<li>Read in original inputs from the repo
<li>Duplicate the input, to ensure there is a reference to comapre to
<li>Check the sequence to see if it contains any non-natural amino acids, mutate back if I do
<li>Check all caps to see if any need to be charged or not. Here teh charged caps were kept, as they look like they interact and neautralise eachother.
<li>Open the Protein preparation wizard (Schrodinger 2024-2), 
<li>Run only step 3 - "Preprocess". The only three options activated are:
   <ul type="disc">
    <li>Convert selomethionines to methionines
    <li>Include peptides when capping termini
    </ul>
<li>Export outputs, named "ligand.sdf", "protein.pdb"
<li>Run input through validation script
<li>Scrub files of sensitive information
</ol>
