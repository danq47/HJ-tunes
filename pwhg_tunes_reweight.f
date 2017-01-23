c dummy subroutine to take in y_higgs and ptj1, and give an output weight.
c Later, we will use something like this to link to SC's interpolating function,
c but for now I'm just trying to see where it fits into the code
      function tune_reweight(yh,ptj1,renscfact_r,facscfact_r)
      implicit none
      real * 8 yh, ptj1, tune_reweight,renscfact_r,facscfact_r
      tune_reweight = 1.0
      if(abs(yh).gt.3) then
         tune_reweight = 0.33
      elseif(abs(yh).gt.1.5) then
         tune_reweight = 3.0
      else
         tune_reweight = 1.0
      endif
      end
