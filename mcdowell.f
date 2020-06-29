program mcdowell

c Modification VI: b = b_o + b_._iso + b_sig >= 0
c                  X = X_sig > = 0
c Restrictions: kappa >= 0
c               -1 <= eta <-- 1
c               0 <= w <= 1

c Declare variables
      integer iSeg, iStep, iStepNum
      real A, alpha, aS, aSDot, aStar, aStarDot
      real B, biso
      real bNot, bNotl, bNot2, bNot3, bNot4, bTotal
      real bsig, eta, bisoDot
      real C, CI, C2, C3, C4, casYN, CStar
      real CTE, D, dbNotdT, dCdT, dHStardT, direc
      real E, El, E2, E3, E4, eDif, eDot, eln
      real elnDot, eOld, eors, eorsEnd, eThDot
      real fixTime
      real HStar, HStarl, HStar2, HStar3, HStar4
      real ovrstrss, meche, MStar, n, Riso
      real numSegs, omegStar, Q, R, rate
      real Rg, sDif, sDot, sOld, strain, stress
      real TD, Temp, TempDot, TempMelt, therme
      real theta, thresh, time, timelnc, timeOld
      real w, X, Xbar, mu, kappa, dXbardT
      real XbarSig, XDot
      real second, LHS, RHS

c Initialize variables
      alpha =0.0
      aS = 0.0
      aStar =0.0
      Riso =0.0
      biso =0.0
      X = 0.0
      eln =0.0
      elnDot =0.0
      eThDot =0.0
      eOld =0.0
      sOld =0.0
      strain = 0.0
      meche = 0.0 
      therme = 0.0
      stress = 0.0
      time =0.0
      timelnc = 0.0
      timeOld =0.0

c Read in constants from data file
      call extract(casYN)
      call extract(TD)
      call extract(Temp)
      call extract(CTE)
      call extract(El)
      call extract(E2)
      call extract(E3)
      call extract(E4)
      call extract(R)
      call extract(HStarl)
      call extract(HStar2)
      call extract(HStar3)
      call extract(HStar4)
      call extract(n)
      call extract(A)
      call extract(B)
      call extract(CI)
      call extract(C2)
      call extract(C3)
      call extract(C4)
      call extract(D)
      call extract(Q)
      call extract(TempMelt)
      call extract(Rg)
      call extract(CStar)
      call extract(MStar)
      call extract(bNotl)
      call extract(bNot2)
      call extract(bNot3)
      call extract(bNot4)
      call extract(eta)
      call extract(kappa)
      call extract(w)
      call extract(mu)
      call extract(thresh)
      call extract(fixTime)
      call extract(numSegs)

c Write out the initial conditions
      write(*,150)time, meche, stress
      do 100 iSeg = 1, int(numSegs)
      
c Read in conditions for each segment
         call extract(TempDot) 
         call extract(eors)
         call extract(rate)
         call extract(eorsEnd)
         
c Determine step times.
         if (abs(rate).It.0.000000001) then
            timelnc = 0.01 * fixTime
            iStep = int(eorsEnd / timelnc)
         elseif (eors.gt.0.0) then
            timelnc = 0.00001 / abs(rate) * fixTime
            iStep = int (abs (log (1.. 0+eorsEnd) - strain) /
     +              abs(rate) / timelnc)
         else
            timelnc = 0.001 / abs(rate) * fixTime
            iStep = int(abs(eorsEnd - stress) /
     +              abs(rate) / timelnc)
         endif

c Run a segment
         do 90 iStepNum = 1, iStep

c Calculate the viscous overstress
            if ((abs(stress - alpha) - R).It.0.0) then
               ovrstrss = 0.0
            else
               ovrstrss = (abs(stress - alpha) - R)
     +                    / (1 - eta)
            endif

c Calculate temperature dependent parameters
            E = El + E2 * tanh(E3 * (Temp - E4))
            HStar = HStarl + HStar2 * tanh(HStar3 * 
     +              (Temp - HStar4))
            dHStardT = TD * HStar2 * HStar3 / (cosh(HStar3
     +                 * (Temp - HStar4))) ** 2.0

c Calculate limiting values of isotropic hardening
            XbarSig = sqrt(2.0/3.0) * kappa * ovrstrss
            Xbar = XbarSig
            
            if (Xbar.gt.0.0000001) then
               XDot = mu * (Xbar - X) * abs(elnDot) + X
     +                / Xbar * dXbardT * TempDot
            else
               XDot = mu * (Xbar - X) * abs(elnDot)
            endif

c Calculate limiting values of backstress
            if(casYN.lt.0.0)then
               bNot = bNotl + bNot2 * tanh(bNot3 *
     +                (Temp - bNot4))

               dbNotdT = TD * bNot2 * bNot3 / (cosh(bNot3
     +                   * (Temp - bNot4))) ** 2.0
            else
               bNot = bNotl + bNot2 * sinh(bNot3 *
     +                (Temp - bNot4))
               dbNotdT = TD * bNot2 * bNot:3 * cosh(bNot3
     +                   * (Temp - bNot4))
            endif
            bsig = sqrt(2.0/3.0) * eta * ovrstrss
            bisoDot = w * XDot
            biso = biso + bisoDot * (timelnc)
            bTotal = bNot + biso + bsig

c Calculate the short range backstress coefficient
            C = CI + C2 * tanh(C3 * (Temp - C4))
            dCdT = TD * C2 * C3 / (cosh(C3 *
     +             (Temp - C4))) ** 2.0

c Calculate theta
            if (Temp.ge.TempMelt / 2.0) then
               theta = exp(-Q / Rg / Temp)
            else
               theta = exp(-2.0 * Q / Rg / TempMelt *
     +                 (log(TempMelt / 2.0 / Temp) + 1.0))
            endif

c Calculate direction
            if(stress.ge.alpha) then
               direc = 1.0
            else
               direc = -1.0
            endif

c Calculate thermal recovery coefficient
c on long-range backstress
            omegStar = CStar * (2.0 / 3.0) ** (MStar / 2.0)
     +                 * aStar ** MStar

c Calculate rates
            if (abs(ovrstrss).It.thresh) then
               elnDot =0.0
               aSDot =0.0
               aStarDot =0.0
               XDot =0.0
            else
            elnDot = A * (ovrstrss / D) ** n * exp(B *
     +               (ovrstrss / D) ** (n + 1.0)) * theta *
     +               direc
            aStarDot = HStar * abs(elnDot) * sqrt(3.0 /
     +                 2.0) * direc + aStar / HStar * dHStardT
     +                 * TempDot -- HStar * omegStar * theta * 
     +                 aStar
            aSDot = C * (sqrt(3.0 / 2.0) * bTotal *
     +              direc - aS) * abs(elnDot) + (aS / C *
     +              dCdT + aS / bTotal * dbNotdT) * TempDot
            endif

            eThDot = CTE * TempDot
            RisoDot = sqrt(3.0/2.0) * (1-w) * XDot

            if (eors.gt.0.0) then
               if (eorsEnd.gt.strain) then
                  eDot = abs(rate) - eThDot
                  sDot = E * (eDot: - elnDot)
               else
                  eDot = -abs(rate) - eThDot
                  sDot = E * (eDot - elnDot)
               endif
            else
               if (eorsEnd.gt.stress) then
                  sDot = abs(rate)
                  eDot = sDot / E + elnDot + eThDot
               else
                  sDot = -abs(rate)
                  eDot = sDot / E + elnDot + eThDot
               endif
            endif

c Calculate new values from rates
            time = time + timelnc
            strain = strain + ((eDot + eThDot) * timelnc)
            eln = eln + (elnDot * timelnc)
            stress = stress + (sDot * timelnc)
            aS = aS + (aSDot * timelnc)
            aStar = aStar + (aStarDot * timelnc)
            alpha = aS + aStar
            X = X + XDot
            R = R + (RisoDot * timelnc)
            Temp = Temp + (TempDot * timelnc)
            therme = therme + (eThDot * timelnc)
            meche = strain - therme

c Check the 2nd Law
            LHS = stress * abs(elnDot) * direc
            if (Xbar.ge.O.00001) then
               RHS = aS * abs(elnDot) * (direc - aS/b) +
     +               aStar * (abs(elnDot) * direc - omegStar
     +               * theta * aStar) + X * abs(elnDot) *
     +               (1 - X/Xbar)
            else
               RHS = aS * abs(elnDot) * (direc - aS/b) +
     +               aStar * (abs(elnDot) * direc - omegStar 
     +               * theta * aStar)
            endif
            second = LHS - RHS

            if (abs(rate).It.0.000000001) then
               if (time-timeOld.gt.5.0) then
                  timeOld = time
                  write(*,200)time, meche, stress,
     +                        second, btotal, Xbar
               endif
            endif

            eDif = abs(eOld - meche)
            sDif = abs(sOld - stress)
c            if (eDif.ge.0.001.or.sDif.ge.0.3) then
            if (eDif.ge.0.00005.or.sDif.gt.0.5) then
               write(*,200)time, meche, stress,
     +                     second, btotal, Xbar
               eOld = meche
               sOld = stress
            endif
90       continue
            write(*,200)time, meche, stress,
     +                  second, btotal, Xbar
100   continue
150   format (fl0.3, fl0.5, fl2.2)
200   format (fl0.3, fl0.5, fl2.2, f12.2, fl2.2, fl2.2
400   format (fl2.2, fl2.2, fl2.2, fl2.2)
      end
c-----------------------------------
      subroutine extract (final)
      integer i, getlen, decim, begin, idirec
      real facten, final, direc
      character *40 indata, newdata
      character *1 blank, period, negative
      parameter (blank = ' ')
      arameter (period = ' . ' )
      parameter (negative = '-' )
      read(*,310) indata
      direc = 1.0
      begin = 0
      decim = 0
      getlen = 0
      idirec = 0
      do 400 i = len(indata), 1, -1
         if(indata(i:i).ne.blank.and.getlen.eq.0)then 
            getlen = i
         endif
         if(indata(i:i).eq.period.and.decim.eq.0)then
            decim = i
         endif
         if(indata(i:i) .eq.negative.and.getlen.ne.0
     +      .and.begin.eq.0.and.idirec.eq.0)then
            begin = i+1
            idirec = i
            direc = -1.0
         endif
         if(indata(i:i).eq.blank.and.getlen.ne.0
     +      .and.begin.eq.0)then
            begin = i+1
            idirec = i
         endif
400   continue

      facten = 0.1
      do 420 i=l, decim-begin
         facten = facten * 10.0
420   continue

      newdata(1:getlen-begin) =
     +     indata(begin:decim-1)//indata(decim+1:getlen)
      final =0.0
      do 410 i=l, getlen-begin
         final = final + (real(ichar(newdata(i:i)))
     +           -48.0)*facten
         facten = facten / 10.0
410   continue

      final = final * direc
      
310   format (a)

return
end 