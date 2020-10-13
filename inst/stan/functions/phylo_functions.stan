  real terms_loglik_lp(real Ne_I_scalar_phy, int M, int N_phy, int N2, vector intervalLength_phy, vector ltt_terms_phy, int[] ne_bin_phy, int[] nco_phy, matrix infectiousness, matrix Rt) {
    vector[N2] ne_binned_phy; // N2 days' worth of Ne values - one bin per day
    real ne_phy; // current value for Ne
    real sterms_phy; // current value for s terms
    real coterms_phy; // current value for co terms
    real lprob_phy; // log prob accumulator
    lprob_phy = 0.0;
    if(N_phy > 0){
      for(j in 1:M){ // eventually loop over M groups. For now M=1.
        for(i in 1:N2){ // write ne as a function of some scalar, I, Rt, and the recovery rate
          ne_binned_phy[i] =  Ne_I_scalar_phy * infectiousness[i,j] / (2.0 * Rt[i,j] * 365.0/6.5 ) ;
        }
      }
      // add up likelihood over all terms
      for(i in 1:N_phy){
        ne_phy = ne_binned_phy[ne_bin_phy[i]]; // ne value corresponding to position i of coalescent terms
        sterms_phy = 0;
        if(ne_phy > 0){
          sterms_phy = -intervalLength_phy[i] * ltt_terms_phy[i]/ne_phy;
          coterms_phy = 0;
          if(ltt_terms_phy[i]>0){
            coterms_phy = nco_phy[i] * (log(ltt_terms_phy[i]) - log(ne_phy));
          }
          lprob_phy = lprob_phy + coterms_phy + sterms_phy;
        }
      }
    }
    return lprob_phy;
  }
