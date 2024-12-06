      PROGRAM EMEGENTGRAVITYMASS

      IMPLICIT DOUBLE PRECISION (A-H,K,M,O-Z) 
      INTEGER I   
      DIMENSION TZERO(12),RT(12),AAKSEN(12),B(12),CAKSEN(12),
     .          TMINTZERO(12),RCOOL(12),ACOOL(12),MNZERO(12),
     .          RC(12),RS(12),ALPHAAKSEN(12),BETA(12),EPSILONN(12),
     .          MNZEROAKSEN(12),RAKSENC(12),BETAAKSEN(12),GAMMA(12),
     .          R(12),RMIN(12),RMAX(12),CMK(12),MNEWT(12)
C--------------------------------------------------------
C result save in file X.dat 
C--------------------
        WRITE(*,*)'-------------------------------------------------'
        WRITE(*,*)'Program to Calculate Mass_Emergent Gravity'
        WRITE(*,*)'-------------------------------------------------'
        WRITE(*,*)' '
        
      
        DO i = 1,12
       !1 cm = 3.2407793D-22 kpc
        !1 cm^(-3) = (3.2407793)**D-66 kpc^(-3)
        !10D-1*1 cm^(-3) = (3.2407793)**(-3)*D65 -> MNZEROAKSEN(CMK)
        !10D-3*1 cm^(-3) = (3.2407793)**(-3)*D63 -> MNZERO(CML)
      

!Code modifikasi gravitasi BHG
        
!Buka file input
      Open (unit=10,file='Data input mass 5Feri rr.dat', status='old')
      
      Read(10,*)TZERO(i),RT(i),AAKSEN(i),B(i),CAKSEN(i),TMINTZERO(i),
     .          RCOOL(i),ACOOL(i),MNZERO(i),RC(i),RS(i),ALPHAAKSEN(i),
     .          BETA(i),EPSILONN(i),MNZEROAKSEN(i),RAKSENC(i),
     .          BETAAKSEN(i),GAMMA(i),R(i),RMIN(i),RMAX(i),CMK(i),
     .          MNEWT(i)
!Seluruh variabel-variabel di persamaan dibuat seperti matriks (dengan label indeks komponen)

      
!Buka file output
        Open (unit=11,file='MASS_emergentGRAV-BHG.dat',action='write')

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
        
!Dari paper, nilai Upsilon, dimensionless, yaitu rentang -0.565 kurang dari sama dengan Upsilon kurang dari sama dengan 0.234 dengan statistik 1sigma (lebih ketat ketelitiannya)
!Sehingga Upsilon adalah nilai rata-rata dari rentang tersebut, (intervalnya 0.799, sehingga mediannya 0.799/2 = 0.3995) 
!Yaitu (-0.565+0.3995) = -0.1655, Upsilon dimensionless dari paper C.Wibisono.
        Upsilon = 1.655D71
!Msun dalam satuan KeV, Msun = 1.11*10^60 MeV = 1.11*10^63 keV
        Msun = 1.11D63 !keV
        MSUNBAR = 10D14*Msun !keV

!Pada BHG terdapat faktor Gn, Gn = G/(1+5*Upsilon). 
        Mstd = (-3.7D-4*T*R(i))*(DLNRHODLNR+DLNTDLNR)
        MMIN = (-3.7D-4*T*RMIN(i))*(DLNRHODLNR+DLNTDLNR)
        MMAX = (-3.7D-4*T*RMAX(i))*(DLNRHODLNR+DLNTDLNR)
        
!Persamaan massa klaster BHG, dengan pendekatan perturbasi M = f + bf^2       
       M = Mstd-(Upsilon*(R(i)**(2)))/(4.0D0*MSUNBAR)*(Mstd**(2.0D0))
       M2= MMIN-(Upsilon*(RMIN(i)**(2)))/(4.0D0*MSUNBAR)*(Mstd**(2.0D0))
       M3= MMAX-(Upsilon*(RMAX(i)**(2)))/(4.0D0*MSUNBAR)*(Mstd**(2.0D0))
       PERSENTASI=(M-Mstd)*100/Mstd
        write(11,*) 'M BHG(dalam 10**14 Msun) = ',M,M2,M3,R(i),RMIN(i),
     &     RMAX(i),PERSENTASI        
        
!Hasilnya aman (tanpa orde tambahan) dan sama dengan code BHG tanpa pendekatan perturbasi, tapi boleh dicek lagi pak 
!Barangkali ada kesalahan

      END DO
      WRITE(*,*)' '
      WRITE(*,*)'Done'
      close(11)
       
      STOP
      
      END PROGRAM
