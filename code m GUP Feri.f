      PROGRAM EMEGENTGRAVITYMASS

      IMPLICIT DOUBLE PRECISION (A-H,K,L,M,O-Z) 
      INTEGER I   
C      Real :: LL, LP
      DIMENSION TZERO(12),RT(12),AAKSEN(12),B(12),CAKSEN(12),
     .          TMINTZERO(12),RCOOL(12),ACOOL(12),MNZERO(12),
     .          RC(12),RS(12),ALPHAAKSEN(12),BETA(12),EPSILONN(12),
     .          MNZEROAKSEN(12),RAKSENC(12),BETAAKSEN(12),GAMMA(12),
     .          R(12),RMIN(12),RMAX(12),CMK(12)
C--------------------------------------------------------
C result save in file X.dat 
C--------------------
        WRITE(*,*)'-------------------------------------------------'
        WRITE(*,*)'Program to Calculate Mass_Emergent Gravity'
        WRITE(*,*)'-------------------------------------------------'
        WRITE(*,*)' '
        
      
        DO i = 1,12
       !1 cm = 3.2407793D-22 kpc
        !1 cm^(-3) = (3.2407793)**(-3)*D66
        !10D-1*1 cm^(-3) = (3.2407793)**(-3)*D65 -> MNZEROAKSEN(CMK)
        !10D-3*1 cm^(-3) = (3.2407793)**(-3)*D63 -> MNZERO(CML)
        
        !CMK = (3.2407793)**(-3)*1.0D65
        !CML = (3.2407793)**(-3)*1.0D63
                  
!Code mass GUP
        
!Buka file input
      Open (unit=10,file='Data input mass 5Feri rr.dat', status='old')
      
      Read(10,*)TZERO(i),RT(i),AAKSEN(i),B(i),CAKSEN(i),TMINTZERO(i),
     .          RCOOL(i),ACOOL(i),MNZERO(i),RC(i),RS(i),ALPHAAKSEN(i),
     .          BETA(i),EPSILONN(i),MNZEROAKSEN(i),RAKSENC(i),
     .          BETAAKSEN(i),GAMMA(i),R(i),RMIN(i),RMAX(i),CMK(i)
!Seluruh variabel-variabel di persamaan dibuat seperti matriks (dengan label indeks komponen)

      
!Buka file output
        Open (unit=11,file='MASS_emergentGRAV-GUP22.dat',action='write')

!Definisi untuk persamaan H
        SRRC = (1.0D0+((R(i)**2)/(RC(i)**2)))
!Sudah di Write(*,*)'SRRC = ',SRRC, hasilnya benar sesuai data paper
     
!Definisi untuk persamaan H
        SRRCA = (1.0D0+(R(i)**2/RAKSENC(i)**2))
!Sudah di Write(*,*)'SRRCA = ',SRRCA, hasilnya benar sesuai data paper
        
!Definisi untuk persamaan H
        RRC = (R(i)/RC(i))
!Sudah di write(*,*)'RRC=',RRC, hasilnya benar sesuai data paper
        
!Definisi untuk persamaan H
        SRGRS = (1.0D0+R(i)**GAMMA(i)*RS(i)**(-GAMMA(i)))
!Sudah di write(*,*)'SRGRS =',SRGRS, hasilnya benar sesuai data paper
        
!Definisi untuk persamaan H
        EPG = (-EPSILONN(i)/GAMMA(i))
!Sudah di write(*,*)'EPG=',EPG, hasilnya benar sesuai data paper
        
!Definisi untuk persamaan H
        PANG = (5.0D-1*ALPHAAKSEN(i)-3.0D0*BETA(i))  
!Sudah di Write(*,*)'PANG =',PANG, hasilnya benar sesuai data paper 

        
!Persamaan suku pertama H
        H1 = (MNZEROAKSEN(i)**2)*SRRCA**(-3.0D0*BETAAKSEN(i))
!Sudah di Write(*,*)'H1 =',H1, hasilnya benar sesuai data paper

!Persamaan suku kedua H
       SRE = (SRGRS**EPG)
      H2 =(MNZERO(i)**2)*(SRRC**PANG)*RRC**(-ALPHAAKSEN(i))*SRE

!Persamaan H
        H = 2.0D0*(H1+H2)

!Definisi untuk persamaan dlnT/dlnr
        RRCAC = (R(i)/RCOOL(i))**ACOOL(i)
!Sudah di write(*,*)'RRCAC = ',RRCAC, hasilnya benar sesuai data paper

!Definisi untuk persamaan dlnT/dlnr
        RRTB = (R(i)/RT(i))**B(i)
!Sudah di Write(*,*)'RRTB =',RRTB, hasilnya benar sesuai data paper
        
!Definisi untuk persamaan dlnT/dlnr
        PYT = (1.0D0+RRCAC)*(TZERO(i)*TMINTZERO(i)+TZERO(i)*RRCAC)
        TMN = (TZERO(i)-TZERO(i)*TMINTZERO(i))/PYT

!Persamaan suku pertama dlnT/dlnr
        DLNTDLNR1 = ACOOL(i)*RRCAC*TMN 

!Persamaan suku kedua dlnT/dlnr = AAKSEN

!Persamaan suku ketiga dlnT/dlnr
        DLNTDLNR3 = (CAKSEN(i)*RRTB)/(1.0D0+RRTB)
        
!Persamaan dlnT/dlnr
        DLNTDLNR = DLNTDLNR1-AAKSEN(i)-DLNTDLNR3
!Sudah di write DLNTDLNR Sesuai
       
!Definisi persamaan dlnrho/dlnr
        MNA = (-MNZERO(i)**2)*ALPHAAKSEN(i)

!Definisi persamaan dlnrho/dlnr
        MNB = MNA*(SRRC**PANG)*(RRC**(-ALPHAAKSEN(i)))*(SRGRS**EPG)

!Definisi persamaan dlnrho/dlnr
        MNC = (-1.0D0+PANG)
        
!Definisi persamaan dlnrho/dlnr
	  MND = (MNZERO(i)**2)*(ALPHAAKSEN(i)-6.0D0*BETA(i))*(SRRC**MNC)
        
!Definisi persamaan dlnrho/dlnr
 	  MNE = RRC**(-ALPHAAKSEN(i)+2.0D0)*(SRGRS**EPG)

!Persamaan suku pertama dlnrho/dlnr
        DLNRHODLNR1 = MNB+MND*MNE
!Sudah sesuai persamaannya

!Definisi persamaan suku kedua dlnrho/dlnr
       RRCA = (R(i)/RAKSENC(i))**(2)
      MNF=6.0D0*(MNZEROAKSEN(i)**2)*BETAAKSEN(i)*RRCA
        MNG = SRRCA**(-1.0D0-3.0D0*BETAAKSEN(i))
!Sudah sesuai persamaannya

!Persamaan suku kedua dlnrho/dlnr
	  DLNRHODLNR2 = MNF*MNG
!Sudah sesuai persamaannya

!Definisi persamaan suku ketiga dlnrho/dlnr
	  MNH = (MNZERO(i)**2)*EPSILONN(i)*(RRC**(-ALPHAAKSEN(i)))
       MNI = (R(i)/RS(i))**(GAMMA(i))*(SRRC**PANG)*(SRGRS**(-1.0D0+EPG))
!Sudah sesuai persamaannya
        
!Persamaan suku ketiga dlnrho/dlnr
		DLNRHODLNR3 = MNH*MNI
!Sudah sesuai persamannya
        
!Persamaan dlnrho/dlnr
	  DLNRHODLNR = (1.0D0/H)*(DLNRHODLNR1-DLNRHODLNR2-DLNRHODLNR3)
!Sudah di write(*,*)'DLNRHODLNR =', DLNRHODLNR, hasil sesuai   

!Definisi persamaan temperature profile
        DDD =(R(i)/RT(i))**(-AAKSEN(i))
!Sudah di Write(*,*)'DDD =',DDD , koreksi tadinya alphaaksen harusnya aaksen.
        
        DDDD = (1.0D0+RRTB)**(CAKSEN(i)/B(i))
!Sudah di write(*,*)'DDDD =',DDDD, hasilnya benar sesuai data paper

        THETAOUT = DDD/DDDD
!Sudah di write(*,*)'Thetaout=',THETAOUT, hasilnya benar sesuai data paper
        
!Definisi persamaan temperature profile
        THETAIN = (RRCAC+TMINTZERO(i))/(RRCAC+1.0D0)
!Sudah di write(*,*)'Thetain = ',THETAIN, hasilnya benar sesuai data paper
        
!Persamaan temperature profile        
        T = TZERO(i)*THETAIN*THETAOUT
!Sudah di write(*,*) 'T =',T, hasilnya benar sesuai data paper
  
!Definisi untuk suku persamaan massa standard
		Mstd = -3.7D-4*T*R(i)*(DLNRHODLNR+DLNTDLNR)
        MMIN = -3.7D-4*T*RMIN(i)*(DLNRHODLNR+DLNTDLNR)
        MMAX = -3.7D-4*T*RMAX(i)*(DLNRHODLNR+DLNTDLNR)
!Sudah di write(*,*)'Mstd=', Mstd, hasilnya benar

		Pi = 3.14

!LP = Panjang Planck yaitu 1.616x10^-35 m, berhubung R(i) nya dalam kpc, maka LP harus diubah ke kpc
!1 m = 3.2407793D-20 kpc
!LP = 1.616x10^-35 m = 1.616*10^-35*(3.2407793D-20) kpc = 5.2370993*D-55 kpc         
        LP = 5.2370993D-55

!Nilai BETANOL didapat dari paper Idrus,(WD and GUP) yaitu BETANOL = -1.656x10^43, ini dimensionless  
		BETANOL = -1.656D43
       
!Definisi untuk suku persamaan massa GUP
		!G = 6.67D-11 !N.m^2/kg^2 
!Satuan dalam SI, sehingga variabelnya ikut menyesuaikan dalam SI, kecuali LP harus diubah ke kpc, karena r dalam kpc
  
   		!Hbar = 1.0545718D-34 !kg.m^2/s
!Msun dalam satuan KeV, Msun = 1.11*10^60 MeV = 1.11*10^63 KeV
        !Msun = 1.99D30 !kg
        !MSUNBAR = Msun/10D14 !kg
        !c = 3.0D8 !m/s
        
        
     	!BETAGUP = BETANOL*(LP**(2))/(Hbar**(2))
!Ketika di write BETAGUP = Nol pak, karena pangkat LP dikuadratkan yang terlalu kecil, anehnya kalau LP dihilangkan, malah keluar hasilnya pak, tidak nol.

!Inputan persamaannya menurut saya sudah sesuai dengan hasil analitik.

!Definisi suku koreksi pada GUP
		BETALP = (BETANOL*(LP**(2)))/(16.0D0*Pi**(2)*R(i)**(2))
        BETALP1=(BETANOL*(LP**(2)))/(16.0D0*Pi**(2)*RMIN(i)**(2))
        BETALP2=(BETANOL*(LP**(2)))/(16.0D0*Pi**(2)*RMAX(i)**(2))
!Sudah di write(*,*)'Betalp=', BETALP,BETALP1,BETALP2, ordenya sangat kecil sekali, sehingga mendekati nol.

        X = (Pi*(R(i)**(2)))/(LP**(2))
        X1 = (Pi*(RMIN(i)**(2)))/(LP**(2))
        X2 = (Pi*(RMAX(i)**(2)))/(LP**(2))
        LL = Log(X)
        LL1 = Log(X1)
        LL2 = Log(X2)
        !Sudah di write(*,*)'LL=',LL,LL1,LL2

        !BETAC = (BETAGUP*G*Hbar*(MSUNBAR**(2)))/(Pi*(R(i)**(2))*c)
        
!Cek hitung ln
        !X = 2
        !LL = Log(X)  !Log = ln
        !Write(*,*)'LL = ', LL
        
!Hasil persamaan massa GUP, diperoleh dengan pendekatan perturbasi M = f+b*f^2 dan ekspansi binomial,
!Nilainya yaitu M = (a-ad)+(-ca(a-ad)^2) = a[1-d-c(a-ad)^2] 
!a = Mstd; d = BETAHA*LL; c = BETAC	      
!Persamaan massa klaster GUP       
        M  = Mstd*((1+BETALP*LL)**(-1.0D0))
        M2 = MMIN*((1+BETALP1*LL1)**(-1.0D0))
        M3 = MMAX*((1+BETALP2*LL2)**(-1.0D0)) 
        write(11,*) 'Mgup(dalam 10**14 Msun) = ', M,M2,M3,R(i),RMIN(i),
     &     RMAX(i)       
!Hasilnya sama dengan Mstandard, karena kontribusi GUP di skala galaksi tidak berpengaruh
        
!Suku koreksi GUP, KI=BETAHA*LL-BETAC*(Mstd-Mstd*BETAHA*LL)**(2)
!Sudah di write(*,*)'KI=',KI, suku kontribusi GUP ini bernilai orde sangat kecil yaitu pangkat di rentang -24 sampai -27,
!Maka hasil massa GUP akan kembali dengan mass standard
        
      END DO
      WRITE(*,*)' '
      WRITE(*,*)'Done'
      close(11)
       
      STOP
      
      END PROGRAM
