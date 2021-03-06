  
       subroutine precip_proc()

	implicit none
	include 'common.inc'
	include 'params.inc'

        integer i,j,k
        real autor, autos, accrr, accris, accrcs, accrig, accrcg
        real evap, dq, omn, omp, omg , qci0
        real pows1, pows2, powg1, powg2, powr1, powr2, tmp
	real qii, qcc, qrr, qss, qgg
        real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
	real f0(nzm),df0(nzm)

        powr1 = (3 + b_rain) / 4.
        powr2 = (5 + b_rain) / 8.
        pows1 = (3 + b_snow) / 4.
        pows2 = (5 + b_snow) / 8.
        powg1 = (3 + b_grau) / 4.
        powg2 = (5 + b_grau) / 8.
      
        if(OCEAN) qci0 = qci0o 
        if(LAND) qci0 = qci0l 
 
	do k=1,nzm
	 do j=1,ny
	  do i=1,nx	  
	  
c-------     Autoconversion/accretion 

	   if(qn(i,j,k)+qp(i,j,k).gt.0.) then

                 omn = omegan(i,j,k)
                 omp = omegap(i,j,k)
                 omg = omegag(i,j,k)

                 qrr = qp(i,j,k) * omp
                 qss = qp(i,j,k) * (1.-omp)*(1.-omg)
                 qgg = qp(i,j,k) * (1.-omp)*omg


		 if(qn(i,j,k).gt.0.) then
     
	           qcc = qn(i,j,k) * omn
	           qii = qn(i,j,k) * (1.-omn)

      	           if(qcc .gt. qcw0) then
      	            autor = alphaelq
!      	            autor = 500.*qcc*qcc   ! Berry Autoconversion (Nc=200)
      	           else
      	            autor = 0.
      	           endif 
      	           if(qii .gt. qci0) then
      	            autos = betaelq*coefice(k)
      	           else
      	            autos = 0.
      	           endif 

     	           accrr = accrrc(k) * qrr ** powr1
     	           tmp = qss ** pows1
     	           accrcs = accrsc(k) * tmp
     	           accris = accrsi(k) * tmp
     	           tmp = qgg ** powg1
     	           accrcg = accrgc(k) * tmp
     	           accrig = accrgi(k) * tmp
                   qcc = (qcc+dtn*autor*qcw0)/
     &			(1.+dtn*(accrr+accrcs+accrcg+autor))
                   qii = (qii+dtn*autos*qci0)/
     &			(1.+dtn*(accris+accrig+autos))
	           dq = dtn *(accrr*qcc + autor*(qcc-qcw0)+
     &		      (accris+accrig)*qii + (accrcs+accrcg)*qcc +
     &                 autos*(qii-qci0))
		   dq = min(dq,qn(i,j,k))
		   qp(i,j,k) = qp(i,j,k) + dq
	           q(i,j,k) = q(i,j,k) - dq
	           qn(i,j,k) = qn(i,j,k) - dq
	           qpsrc(k) = qpsrc(k) + dq


	         elseif(qp(i,j,k).gt.qp_threshold
     &                 .and.qn(i,j,k).eq.0.) then
	           dq = dtn * (evapr1(k) * sqrt(qrr) + 
     &                		evapr2(k) * qrr**powr2 +
     &				evaps1(k) * sqrt(qss) +
     &				evaps2(k) * qss**pows2 +
     &				evapg1(k) * sqrt(qgg) +
     &				evapg2(k) * qgg**powg2)*
     &			 	(q(i,j,k) /qsat(i,j,k)-1.) 
                   dq = max(-0.5*qp(i,j,k),dq) 
                   qp(i,j,k) = qp(i,j,k) + dq
                   q(i,j,k) = q(i,j,k) - dq
	           qpevp(k) = qpevp(k) + dq

		 else
	
	           q(i,j,k) = q(i,j,k) + qp(i,j,k)
	           qp(i,j,k) = 0.
	           qpevp(k) = qpevp(k) - qp(i,j,k)

                 endif

	    endif

      	    qp(i,j,k)=max(0.,qp(i,j,k))

	  end do
	 enddo
	enddo
    


      return
      end

