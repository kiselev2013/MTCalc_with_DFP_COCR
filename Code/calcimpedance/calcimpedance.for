      implicit real*4 (a-h,o-z)      
	parameter(idim=20000) 
	dimension B(4,4), Q(4,4), F(4), U(4), F1(4), OmegaM(idim),
     *Omegat(idim),Ud(4)
	dimension Xkoor(idim),Ykoor(idim),Zkoor(idim)
	pi=3.14159265358979323846264
	omu0=4d0*pi*1d-7
	open(1,file='frec')
	read(1,*) KolOmega
	do i=1,KolOmega
	read(1,*) Omegat(KolOmega-i+1)
	OmegaM(KolOmega-i+1)=2d0*pi*Omegat(KolOmega-i+1)
	enddo
	close(1)

	open(2,file='sqrt(T)')
	write(2,*) KolOmega,'   sqrt(T)     nu'
	do i=KolOmega,1,-1
	write(2,*) sqrt(1./Omegat(i)),Omegat(i)
	
	enddo
	close(2)
	

	open(1,file='pointres')
	read(1,*) kol
	open(2,file='ex_sx_all_s0')
	open(3,file='ex_cx_all_s0')
	open(4,file='ey_sx_all_s0')
	open(5,file='ey_cx_all_s0')
	open(6,file='hx_sx_all_s0')
	open(7,file='hx_cx_all_s0')
	open(8,file='hy_sx_all_s0')
	open(9,file='hy_cx_all_s0')
	open(10,file='hz_sx_all_s0')
	open(11,file='hz_cx_all_s0')
	open(12,file='ex_sy_all_s0')
	open(13,file='ex_cy_all_s0')
	open(14,file='ey_sy_all_s0')
	open(15,file='ey_cy_all_s0')
	open(16,file='hx_sy_all_s0')
	open(17,file='hx_cy_all_s0')
	open(18,file='hy_sy_all_s0')
	open(19,file='hy_cy_all_s0')
	open(58,file='hz_sy_all_s0')
	open(59,file='hz_cy_all_s0')

	call openres()


	do m=1,kol
	read(1,*) Xkoor(m),Ykoor(m),Zkoor(m)

	do jjj=1,1

	do j1=2,19
	call readgl(j1)
	enddo
	call readgl(58)
	call readgl(59)

	do j2=1,KolOmega

	omega=OmegaM(j2)
	read(2,*) x,ex_s0
	read(3,*) x,ex_c0
	read(4,*) x,ey_s0
	read(5,*) x,ey_c0
	read(6,*) x,hx_s0
	read(7,*) x,hx_c0
	read(8,*) x,hy_s0
	read(9,*) x,hy_c0
	read(10,*) x,hz_s0
	read(11,*) x,hz_c0
	read(12,*) x,ex_s1
	read(13,*) x,ex_c1
	read(14,*) x,ey_s1
	read(15,*) x,ey_c1
	read(16,*) x,hx_s1
	read(17,*) x,hx_c1
	read(18,*) x,hy_s1
	read(19,*) x,hy_c1
	read(58,*) x,hz_s1
	read(59,*) x,hz_c1


	if (jjj.eq.1) then

	B(1,1)=hx_s0
	B(1,2)=-hx_c0
	B(1,3)=hy_s0
	B(1,4)=-hy_c0
	B(2,1)=hx_c0
	B(2,2)=hx_s0
	B(2,3)=hy_c0
	B(2,4)=hy_s0
	B(3,1)=hx_s1
	B(3,2)=-hx_c1
	B(3,3)=hy_s1
	B(3,4)=-hy_c1
	B(4,1)=hx_c1
	B(4,2)=hx_s1
	B(4,3)=hy_c1
	B(4,4)=hy_s1

	call  makeQ (B,Q)


	F(1)=ex_s0
	F(2)=ex_c0
	F(3)=ex_s1
	F(4)=ex_c1
	do i=1,4
	U(i)=0d0
	do j=1,4
	U(i)=U(i)+Q(i,j)*F(j)
	enddo
	enddo



	Zxx=sqrt(U(1)*U(1)+U(2)*U(2))
	if (U(1).gt.0..and.U(2).gt.0.) Fxx=atan(U(2)/U(1))
	if (U(1).lt.0..and.U(2).gt.0.) Fxx=atan(U(2)/U(1))+pi
	if (U(1).lt.0..and.U(2).lt.0.) Fxx=atan(U(2)/U(1))-pi
	if (U(1).gt.0..and.U(2).lt.0.) Fxx=atan(U(2)/U(1))

	if (U(1).gt.0..and.abs(U(2)).lt.1e-15) Fxx=0.
	if (abs(U(1)).lt.1e-15.and.U(2).gt.0.) Fxx=pi/2.
	if (U(1).lt.0..and.abs(U(2)).lt.1e-15) Fxx=pi
	if (abs(U(1)).lt.1e-15.and.U(2).lt.0.) Fxx=-pi/2.

	write(60,rec=KolOmega*m-j2+1) U(1)
	write(61,rec=KolOmega*m-j2+1) U(2)

	Zxy=sqrt(U(3)*U(3)+U(4)*U(4))
	if (U(3).gt.0..and.U(4).gt.0.) Fxy=atan(U(4)/U(3))
	if (U(3).lt.0..and.U(4).gt.0.) Fxy=atan(U(4)/U(3))+pi
	if (U(3).lt.0..and.U(4).lt.0.) Fxy=atan(U(4)/U(3))-pi
	if (U(3).gt.0..and.U(4).lt.0.) Fxy=atan(U(4)/U(3))

	if (U(3).gt.0..and.abs(U(4)).lt.1e-15) Fxy=0.
	if (abs(U(3)).lt.1e-15.and.U(4).gt.0.) Fxy=pi/2.
	if (U(3).lt.0..and.abs(U(4)).lt.1e-15) Fxy=pi
	if (abs(U(3)).lt.1e-15.and.U(4).lt.0.) Fxy=-pi/2.

	write(62,rec=KolOmega*m-j2+1) U(3)
	write(63,rec=KolOmega*m-j2+1) U(4)



	F(1)=ey_s0
	F(2)=ey_c0
	F(3)=ey_s1
	F(4)=ey_c1
	do i=1,4
	Ud(i)=U(i)
	U(i)=0d0
	do j=1,4
	U(i)=U(i)+Q(i,j)*F(j)
	enddo
	enddo


	Zyx=sqrt(U(1)*U(1)+U(2)*U(2))
	if (U(1).gt.0..and.U(2).gt.0.) Fyx=atan(U(2)/U(1))
	if (U(1).lt.0..and.U(2).gt.0.) Fyx=atan(U(2)/U(1))+pi
	if (U(1).lt.0..and.U(2).lt.0.) Fyx=atan(U(2)/U(1))-pi
	if (U(1).gt.0..and.U(2).lt.0.) Fyx=atan(U(2)/U(1))

	if (U(1).gt.0..and.abs(U(2)).lt.1e-15) Fyx=0.
	if (abs(U(1)).lt.1e-15.and.U(2).gt.0.) Fyx=pi/2.
	if (U(1).lt.0..and.abs(U(2)).lt.1e-15) Fyx=pi
	if (abs(U(1)).lt.1e-15.and.U(2).lt.0.) Fyx=-pi/2.

	write(64,rec=KolOmega*m-j2+1) U(1)
	write(65,rec=KolOmega*m-j2+1) U(2)



	Zyy=sqrt(U(3)*U(3)+U(4)*U(4))
	if (U(3).gt.0..and.U(4).gt.0.) Fyy=atan(U(4)/U(3))
	if (U(3).lt.0..and.U(4).gt.0.) Fyy=atan(U(4)/U(3))+pi
	if (U(3).lt.0..and.U(4).lt.0.) Fyy=atan(U(4)/U(3))-pi
	if (U(3).gt.0..and.U(4).lt.0.) Fyy=atan(U(4)/U(3))

	if (U(3).gt.0..and.abs(U(4)).lt.1e-15) Fyy=0.
	if (abs(U(3)).lt.1e-15.and.U(4).gt.0.) Fyy=pi/2.
	if (U(3).lt.0..and.abs(U(4)).lt.1e-15) Fyy=pi
	if (abs(U(3)).lt.1e-15.and.U(4).lt.0.) Fyy=-pi/2.

	write(66,rec=KolOmega*m-j2+1) U(3)
	write(67,rec=KolOmega*m-j2+1) U(4)


	write(20,rec=KolOmega*m-j2+1) Zxx
	write(21,rec=KolOmega*m-j2+1) -Fxx*180/pi
	write(22,rec=KolOmega*m-j2+1) Zxy
	write(23,rec=KolOmega*m-j2+1) -Fxy*180/pi
	write(24,rec=KolOmega*m-j2+1) Zyx
	write(25,rec=KolOmega*m-j2+1) -Fyx*180/pi
	write(26,rec=KolOmega*m-j2+1) Zyy
	write(27,rec=KolOmega*m-j2+1) -Fyy*180/pi


	write(32,rec=KolOmega*m-j2+1) Zxx*Zxx/(omu0*omega)
	write(33,rec=KolOmega*m-j2+1) -2.*Fxx*180/pi
	write(34,rec=KolOmega*m-j2+1) Zxy*Zxy/(omu0*omega)
	write(35,rec=KolOmega*m-j2+1) -2.*Fxy*180/pi
	write(36,rec=KolOmega*m-j2+1) Zyx*Zyx/(omu0*omega)
	write(37,rec=KolOmega*m-j2+1) -2.*Fyx*180/pi
	write(38,rec=KolOmega*m-j2+1) Zyy*Zyy/(omu0*omega)
	write(39,rec=KolOmega*m-j2+1) -2.*Fyy*180/pi




c рассчет эффективных значений
	zeff1d=Ud(1)*U(3)-Ud(2)*U(4)
	zeff1m=Ud(1)*U(4)+Ud(2)*U(3)

	zeff2d=U(1)*Ud(3)-U(2)*Ud(4)
	zeff2m=U(1)*Ud(4)+U(2)*Ud(3)

	zefffd=zeff1d-zeff2d
	zefffm=zeff1m-zeff2m
	
	zmod=sqrt(sqrt(zefffd*zefffd+zefffm*zefffm))
	if (zefffd.gt.0..and.zefffm.gt.0.) zfi=atan(zefffm/zefffd)
	if (zefffd.lt.0..and.zefffm.gt.0.) zfi=atan(zefffm/zefffd)+pi
	if (zefffd.lt.0..and.zefffm.lt.0.) zfi=atan(zefffm/zefffd)-pi
	if (zefffd.gt.0..and.zefffm.lt.0.) zfi=atan(zefffm/zefffd)
	zfi=zfi/2

	write(40,rec=KolOmega*m-j2+1) Zmod*Zmod/(omu0*omega)
	write(41,rec=KolOmega*m-j2+1) -2.*zfi*180/pi
	


c матрица Визе-Паркинсона

	if (m.eq.1108) then
	rrr=1
	endif


	F(1)=hz_s0
	F(2)=hz_c0
	F(3)=hz_s1
	F(4)=hz_c1
	do i=1,4
	U(i)=0d0
	do j=1,4
	U(i)=U(i)+Q(i,j)*F(j)
	enddo
	enddo


	Wzx=sqrt(U(1)*U(1)+U(2)*U(2))
	if (U(1).gt.0..and.U(2).gt.0.) WFzx=atan(U(2)/U(1))
	if (U(1).lt.0..and.U(2).gt.0.) WFzx=atan(U(2)/U(1))+pi
	if (U(1).lt.0..and.U(2).lt.0.) WFzx=atan(U(2)/U(1))-pi
	if (U(1).gt.0..and.U(2).lt.0.) WFzx=atan(U(2)/U(1))

	if (U(1).gt.0..and.abs(U(2)).lt.1e-15) WFzx=0.
	if (abs(U(1)).lt.1e-15.and.U(2).gt.0.) WFzx=pi/2.
	if (U(1).lt.0..and.abs(U(2)).lt.1e-15) WFzx=pi
	if (abs(U(1)).lt.1e-15.and.U(2).lt.0.) WFzx=-pi/2.

	Wzy=sqrt(U(3)*U(3)+U(4)*U(4))
	if (U(3).gt.0..and.U(4).gt.0.) WFzy=atan(U(4)/U(3))
	if (U(3).lt.0..and.U(4).gt.0.) WFzy=atan(U(4)/U(3))+pi
	if (U(3).lt.0..and.U(4).lt.0.) WFzy=atan(U(4)/U(3))-pi
	if (U(3).gt.0..and.U(4).lt.0.) WFzy=atan(U(4)/U(3))

	if (U(3).gt.0..and.abs(U(4)).lt.1e-15) WFzy=0.
	if (abs(U(3)).lt.1e-15.and.U(4).gt.0.) WFzy=pi/2.
	if (U(3).lt.0..and.abs(U(4)).lt.1e-15) WFzy=pi
	if (abs(U(3)).lt.1e-15.and.U(4).lt.0.) WFzy=-pi/2.

	write(68,rec=KolOmega*m-j2+1) U(1)
	write(69,rec=KolOmega*m-j2+1) U(2)
	write(70,rec=KolOmega*m-j2+1) U(3)
	write(71,rec=KolOmega*m-j2+1) U(4)



	write(28,rec=KolOmega*m-j2+1) Wzx
	write(29,rec=KolOmega*m-j2+1) -WFzx*180./pi
	write(30,rec=KolOmega*m-j2+1) Wzy
	write(31,rec=KolOmega*m-j2+1) -WFzy*180./pi

	endif

	enddo

	enddo

	enddo
c
c
111	call closeres()

	end
C
C
C

C
C
      subroutine makeQ (B,Q)
      implicit real*4 (a-h,o-z)      
      dimension B(4,4), Q(4,4), bminor(3,3)
      fdet=0d0
      call determinan(B,fdet)
      do 13 i=1,4
      do 14 j=1,4
      i1=1
      do 11 k1=1,4
      i2=1
      if (k1.ne.i) then
      do 12 k2=1,4
      if (k2.ne.j) then
      bminor(i1,i2)=B(k1,k2)
      i2=i2+1
      end if
   12 continue
      i1=i1+1
      end if
   11 continue
      call minor(bminor,qm)
      Q(j,i)=qm
      s=1
      do 20 m=1,i+j
   20 s=s*(-1)
      Q(j,i)=Q(j,i)*s/fdet
   14 continue
   13 continue
      return
      end
c
      subroutine minor(bm,qm)
      implicit real*4 (a-h,o-z)      
      dimension bm(3,3)
      qm=bm(1,1)*bm(2,2)*bm(3,3)+bm(2,1)*bm(3,2)*bm(1,3)+
     *bm(1,2)*bm(2,3)*bm(3,1)-bm(1,3)*bm(2,2)*bm(3,1)-
     *bm(1,2)*bm(2,1)*bm(3,3)-bm(1,1)*bm(2,3)*bm(3,2)
      return
      end
c
      subroutine determinan(B,fdet)
      implicit real*4 (a-h,o-z)      
      dimension B(4,4),bminor(3,3),Q(4,4)
      fdet=0d0
      i=1
      do 14 j=1,4
      i1=1
      do 11 k1=1,4
      i2=1
      if (k1.ne.i) then
      do 12 k2=1,4
      if (k2.ne.j) then
      bminor(i1,i2)=B(k1,k2)
      i2=i2+1
      end if
   12 continue
      i1=i1+1
      end if
   11 continue
      call minor(bminor,qm)
      Q(i,j)=qm
      s=1
      do 20 m=1,i+j
   20 s=s*(-1)
      Q(i,j)=Q(i,j)*s
      fdet=B(i,j)*Q(i,j)+fdet
   14 continue
      return
      end
c
c
c
	subroutine readgl(j)
      implicit real*4 (a-h,o-z)      
	do k=1,6
	read(j,*)
	enddo
	end
c
c
	subroutine openres()
	open(20,file='Zxx',access='direct',recl=4)
	open(21,file='Fxx',access='direct',recl=4)
	open(22,file='Zxy',access='direct',recl=4)
	open(23,file='Fxy',access='direct',recl=4)
	open(24,file='Zyx',access='direct',recl=4)
	open(25,file='Fyx',access='direct',recl=4)
	open(26,file='Zyy',access='direct',recl=4)
	open(27,file='Fyy',access='direct',recl=4)
	open(28,file='Wzx',access='direct',recl=4)
	open(29,file='WFzx',access='direct',recl=4)
	open(30,file='Wzy',access='direct',recl=4)
	open(31,file='WFzy',access='direct',recl=4)

	open(32,file='Rxx',access='direct',recl=4)
	open(33,file='Fixx',access='direct',recl=4)
	open(34,file='Rxy',access='direct',recl=4)
	open(35,file='Fixy',access='direct',recl=4)
	open(36,file='Ryx',access='direct',recl=4)
	open(37,file='Fiyx',access='direct',recl=4)
	open(38,file='Ryy',access='direct',recl=4)
	open(39,file='Fiyy',access='direct',recl=4)
	open(40,file='Reff',access='direct',recl=4)
	open(41,file='Fieff',access='direct',recl=4)



	open(60,file='ZxxRe',access='direct',recl=4)
	open(61,file='ZxxIm',access='direct',recl=4)
	open(62,file='ZxyRe',access='direct',recl=4)
	open(63,file='ZxyIm',access='direct',recl=4)
	open(64,file='ZyxRe',access='direct',recl=4)
	open(65,file='ZyxIm',access='direct',recl=4)
	open(66,file='ZyyRe',access='direct',recl=4)
	open(67,file='ZyyIm',access='direct',recl=4)
	open(68,file='WzxRe',access='direct',recl=4)
	open(69,file='WzxIm',access='direct',recl=4)
	open(70,file='WzyRe',access='direct',recl=4)
	open(71,file='WzyIm',access='direct',recl=4)



	return
	end

c
	subroutine closeres()
	do j=20,39
	close(j)
	enddo
	end