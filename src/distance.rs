/// calculates the Levenshtein edit distance between two byte slices (kmers)
/// (optimized using a two-row matrix approach to save memory)
pub fn edit_distance(s1: &[u8], s2: &[u8]) -> usize {
    // optimization: ensure s2 is the shorter slice to minimize the memory
    // allocated for the rows.
    // TODO: check that it is correctly
    let (s1, s2) = if s1.len() < s2.len() {
        (s2, s1)
    } else {
        (s1, s2)
    };

    if s2.is_empty() {
        return s1.len();
    }

    // allocate only two rows instead of a full NxM matrix
    // previous_row is initialized with 0, 1, 2, 3...
    let mut previous_row: Vec<usize> = (0..=s2.len()).collect();
    let mut current_row = vec![0; s2.len() + 1];

    for (i, &c1) in s1.iter().enumerate() {
        current_row[0] = i + 1;

        for (j, &c2) in s2.iter().enumerate() {
            let insertions = previous_row[j + 1] + 1;
            let deletions = current_row[j] + 1;

            // branchless substitution cost: evaluates to 1 if different, 0 if same
            let substitutions = previous_row[j] + usize::from(c1 != c2);

            current_row[j + 1] = insertions.min(deletions).min(substitutions);
        }

        // swap the current row into the previous row for the next iteration
        previous_row.copy_from_slice(&current_row);
    }

    previous_row[s2.len()]
}