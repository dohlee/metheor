


struct CpGState {
    pos: i32,
    methylated: bool
}

impl CpGState {
    fn new(pos: i32, c: char) -> CpGState {
        CpGState { pos: pos, methylated: c == 'Z' }
    }
}

impl ToString for CpGState {
    fn to_string(&self) -> String {
        format!("{}, {}", self.pos, self.methylated)
    }
}

pub fn count_z(meth_str: &str) -> i32 {
    (meth_str.matches("z").count() + meth_str.matches("Z").count()) as i32
}

pub fn compute_pairwise_concordance_discordance(cpgs: Vec<(i32, char)>, min_distance: i32, max_distance: i32) -> (i32, i32) {
    let mut anchors: Vec<CpGState> = Vec::new();
    let mut min_anchor_pos = -1;
    let mut n_concordant = 0;
    let mut n_discordant = 0;

    // for (i, c) in meth_str.chars().enumerate().filter(|(_i, c)| (*c == 'z') || (*c == 'Z')) {
    for (i, c) in cpgs.iter(){
        let cpg_state = CpGState::new(*i as i32, *c);

        if min_anchor_pos != -1 {
            while (cpg_state.pos - min_anchor_pos > max_distance) && (anchors.len() > 0) {
                anchors.remove(0);

                if anchors.len() > 0 {
                    min_anchor_pos = anchors[0].pos;
                } else {
                    min_anchor_pos = -1;
                }
            }
        }

        for anchor in anchors.iter() {
            if cpg_state.pos - anchor.pos < min_distance {
                continue;
            }
            
            if anchor.methylated == cpg_state.methylated {
                n_concordant += 1;
            } else {
                n_discordant += 1;
            }
        }

        if min_anchor_pos == -1 { min_anchor_pos = cpg_state.pos; }
        anchors.push(cpg_state);
    }

    // println!("{}, {}", n_concordant, n_discordant);
    (n_concordant, n_discordant)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairwise_concordance_discordance_no_pair() {
        let min_distance = 2;
        let max_distance = 16;
        let tmp_cpgs: Vec<(i32, char)> = Vec::from([(2, 'z'), (100, 'z')]);
        let (n_con, n_dis) = compute_pairwise_concordance_discordance(tmp_cpgs, min_distance, max_distance);

        assert_eq!(n_con, 0);
        assert_eq!(n_dis, 0);
    }

    #[test]
    fn pairwise_concordance_discordance_complex() {
        let min_distance = 2;
        let max_distance = 16;
        let tmp_cpgs: Vec<(i32, char)> = Vec::from([(2, 'z'), (8, 'Z'), (20, 'z'), (32, 'z'), (34, 'z'), (100, 'Z'), (116, 'z')]);
        let (n_con, n_dis) = compute_pairwise_concordance_discordance(tmp_cpgs, min_distance, max_distance);

        assert_eq!(n_con, 3);
        assert_eq!(n_dis, 3);
    }
}
