program bp

C  Define the integers
      integer bigloop, dseg, stepnumr, numline, j

C  Define the real numbers
      real time, stress, stresdot
      real strain, eIN, edot, eINdot
      real Temp, Tempdot
      real SScontrl, SSrate, SSlimit
      real elasMod, n, Al, A2, ZO, Zl, Z2, Z3
      real temper(lOO) , nn(100) , AAl(lOO), ZZO(IOO)
      real iranl(100), mm2(100), EE(100) , ZZ3(100 )
      real ml , m2, rl , r 2
      real steptime, segments
      real dt, ZIS, ZD, ZlSdot, ZDdot
      real dstress, dstrain, sCLD, eOLD
      real sinStres, sinbeta
      real slopeE, slopeZO, slopeZ3, betadot, beta

C  Define initial conditions
      time = 0.0
      stress = 0.0
      strain = 0.0
      eIN =0.0
      ZD = 0.0
      beta =0.0
      betadot =0.0
      eINdot =0.0
      ZlSdot =0.0
      ZDdot =0.0
      SOLD =0.0
      eOLD =0.0

C  Writes the included variables to (*) which
C     is later directed to output file
C
C  Formats written on line 200

      write(*,2 00)strain,stress,ZIS,ZD

C  Extract input from data file
      call extract(Temp)
      call extract(Zl)
      call extract(rl)
      r2 = rl 
      call extract(steptime)
      call extract(segments)

C  Import lookup table data from the lookup file

      call import (numline, temper, EE, rm, AAl, ZZO,
     +             ZZ3, mml, mm2)
      eINdot =0.0
      do 90 bigloop = 1, int(segments)

C  These variables set testing parameters for temperature rate of change,
C  stress/strain control (-1,1), stress/strain rate, and segment limit for
C  stress/strain. We extract them from the data file input.

         call extract(Tempdot)
         call extract(SScontrl)
         call extract(SSrate)
         call extract(SSlimit)

C  Now we normalize dt, the time increment, such that the same number of
steps, dseg,
C  occur in each constant strain rate segment.

         if (abs(SSrate).It.0.000000001) then
            dt = 0.01 * steptime
            dseg = int(SSlimit / dt)
         elseif (SScontrl.gt.0.0) then
            dt = 0.00001 / abs(SSrate) * steptime
            dseg = int(SSlimit / abs(SSrate) / dt)
         else
            dt = 0.001 / abs(SSrate) * steptime
            dseg = int(abs(SSlimit - stress) /
     +             abs(SSrate) / dt)
         endif

C  Now we solve the basic equations for the Bodner-Partom model....

      do 100 stepnumr = 1, dseg

C  The elastic modulus is defined by:
C
C  elasMod = 46000 - 333 * exp(.0112 * (Temp))
C
C  and we use the modulus values solved for at each lookup table
C  data point and then interpolate between these values to solve for
C  any given temperature. This is a more computationally efficient
C  approximation for a parameter that is not a critial determinant
C  of material behavior.
C
C  We now interpolate the temperature dependent constants for our current 
C  temperature

            do 150 j = 1, numline-1
               if (Temp.gt.(temper(j) - 0.0001).and.
     +         temp.It.temper(j+1)) then
                  Al = AAl(j) + (AAKj+1) - AAl
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j)
                  n = nn(j) + (nn(j +1) - nn(j) )
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j)
                  Z0 = ZZO(j) + (ZZ0(j+l) - ZZO
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j)
                  ml = mml(j) + (mml(j+l) - mml
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j)
                  m2 = mm2(j) + (mm2(j+l) - mm2
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j)
                  elasMod = EE(j) + (EE(j+l) - EE(j))*
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j))/
                  Z3 = ZZ3(j) + (ZZ3(j+l) - ZZ3(j))*
     +                 (Temp - temper(j))/
     +                 (temper(j+1) - temper(j))
                  slopeE = (EE(j+l) - EE(j))/
     +                 (temper(j+1) - temper(j))
                  slopeZO = (ZZO(j+1) - ZZ0(j))/
     +                 (temper(j+1) - temper(j))
                  slopeZ3 = (ZZ3(j+l) - ZZ3(j))/
     +                 (temper(j+1) - temper(j))
               endif
150         continue
            A2 = Al
            Z2 = Z0
            if (stepnumr.eq.1.and.
     +         bigloop.eq.1) then
               ZIS = Z0
            endif

C  We now simply define u and v, the sign of the stress and
C  the sign of beta as sinStres and sinbeta, respectively

            if (abs(stress) .It .0.000000001) then
               sinStres = 1.0
            else
               sinStres = stress / abs(stress)
            endif

            if (abs(beta).It.0.000000001) then 
               sinbeta = 1.0
            else
               sinbeta = beta / abs(beta)
            endif

C  Solving the differential BP equations using forward Euler method...

            if (abs(stress).gt.0.000000001) then
               eINdot = 1.1547*10000.0*sinStres*
     +               exp(-0.5*((ZIS+ZD)/
     +               abs(stress))**(2.0*n))
            else
               eINdot = 0
            endif

            ZlSdot = ml*(Zl-ZIS)*stress*eINdot-
     +               A1*Z1*((ZIS-Z2)/Zl)**{rl)+
     +               (Zl-ZIS)/(Z1-Z2)*slopeZ0*Tempdot

            betadot = m2*(Z3*sinStres-beta)*stress*eINdot
     +               -A2*Z1*((abs(beta)/Zl)**r2)*sinbeta
     +               + beta/Z3*slopeZ3*Tempdot

C  Define the experiment as a strain controlled or stress
C  controlled process by reading the "SS" input

            if (SScontrl.gt.0.0) then
               edot = SSrate
               stresdot = elasMod * (edot - eINdot) +
     +             slopeE * Tempdot * (strain - eIN)
            else
               if (SSlimit.gt.stress) then
                  stresdot = abs(SSrate)
                  edot = stresdot / elasMod + eINdot
               else
                  stresdot = -abs(SSrate)
                  edot = stresdot / elasMod + eINdot
               endif
            endif

C  Now we increment our variables by dt*rate

            time = time + dt
            strain = strain + (edot * dt)
            eIN = eIN + (eINdot * dt)
            ZIS = ZIS + (ZlSdot * dt)
            beta = beta + (betadot * dt)
            ZD = beta * sinStres
            stress = stress + (stresdot * dt)
            Temp = Temp + (Tempdot * dt) 

C  Now we select which points are significantly different from the
C  last plotted point so that these points in turn may be plotted. The
C  selected points are written to the output data file.

            dstrain = abs(eOLD - strain)
            dstress = abs(sOLD - stress)
            if (dstrain.ge.0.001.or.dstress.ge.0.2) then
               write(*,200)strain, stress, ZIS, ZD
               eOLD = strain
               sOLD = stress
            endif

100      continue

C  This line plots the very last point for every constant strain segment
C  to ensure "crisp" transition to next segment.

            write(*,200)strain, stress, ZIS,ZD
90    continue
200   format (fl2.8, fl2.4, fl2.4, fl2.4, fl2.4)
      end
C  ----------------------------------
      subroutine extract (final)
      integer i, getlen, decim, begin, idirec
      real facten, final, direc
      character *40 indata, newdata
      character *1 blank, period, negative
      parameter (blank = ' ')
      parameter (period = '.')
      parameter (negative = '- ' )

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
     +     .and.begin.eq.0)then
            begin = i+1
            idirec = i
         endif
400   continue

      facten =0.1
      do 420 i=l, decim-begin
         facten = facten * 10.0
420   continue

      newdata (1 .-getlen-begin) =
     +     indata(begin:decim-1)//indata(decim+1:getlen)

      final =0.0
      do 410 i=l, getlen-begin
         final = final f (real(ichar(newdata(i:i)))
     +            -48.0)*facten
         facten = facten / 10.0
410   continue

      final = final * direc

310   format (a)

      return
      end

C  ----------------------------------
      subroutine import(numline,temp,elasMod,n,A,Z0,Z3,ml,m2)
      integer numline, i
      real temp(lOO), elasMod(lOO), n(100), A(100), Z0(100),
     +         Z3 (100) , ml(100), m2(100)
      open (unit = 1, file = 'lookup', status = 'old')
      read (1, *) numline
      do 100 i = 1, numline
         read (1, 200) temp(i), elasMod(i), n(i), A(i),
     +         Z0(i), Z3 (i), ml(i), m2(i)
100   continue

      close ( Unit=l )
      
200   format (f5.1, 3X, f7.1, 3X, £5.3, 3X, f5.3, 3X,
     +         f5.1, 3X, f5.1, 3X, £5.2, 3X, f5.1)
      return
end