/// Calculates the weighted edit distance between two byte slices (kmers).
/// Allows independent penalization of insertions, deletions, and substitutions.
/// (Optimized using a two-row matrix approach to save memory)
pub fn edit_distance_weighted(
    s1: &[u8],
    s2: &[u8],
    ins_cost: usize,
    del_cost: usize,
    sub_cost: usize,
) -> usize {
    // Optimization: ensure s2 is the shorter slice to minimize memory allocation.
    // CRITICAL FIX: If we swap the sequences, we MUST also swap the insertion
    // and deletion costs to maintain the biological direction of the alignment!
    let (s1, s2, active_ins_cost, active_del_cost) = if s1.len() < s2.len() {
        (s2, s1, del_cost, ins_cost)
    } else {
        (s1, s2, ins_cost, del_cost)
    };

    if s2.is_empty() {
        return s1.len() * active_del_cost;
    }

    // Initialize previous_row with cumulative insertion costs 
    // (Instead of 0, 1, 2... it becomes 0, ins, 2*ins...)
    let mut previous_row: Vec<usize> = (0..=s2.len()).map(|x| x * active_ins_cost).collect();
    let mut current_row = vec![0; s2.len() + 1];

    for (i, &c1) in s1.iter().enumerate() {
        // The first column represents cumulative deletion costs
        current_row[0] = (i + 1) * active_del_cost;

        for (j, &c2) in s2.iter().enumerate() {
            let insertions = previous_row[j + 1] + active_ins_cost;
            let deletions = current_row[j] + active_del_cost;

            // Apply substitution penalty if bases differ, otherwise cost is 0
            let substitutions = previous_row[j] + if c1 != c2 { sub_cost } else { 0 };

            current_row[j + 1] = insertions.min(deletions).min(substitutions);
        }

        // Swap the current row into the previous row for the next iteration
        previous_row.copy_from_slice(&current_row);
    }

    previous_row[s2.len()]
}