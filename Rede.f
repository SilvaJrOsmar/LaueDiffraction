      PROGRAM REDE
C
C-----------------------------------------------------------------------
C     SIMULACAO COMPUTACIONAL DA REDE CRISTALINA A PARTIR
C     DAS MANCHAS DE LAUE
C
C     Calcula n(r), densidade eletronica numa grade 3D
C
C     8: INPUT  Pontos ou manchas de Laue no anteparo 'fig_laue.dat'
C     9: OUTPUT Densidade n(r) para cada x,y,z 'rede.dat'
C
C     Osmar S.Silva Jr. , 2017/09/19
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(300),Y(300)
      real etime
      real elapsed(2)
      real total
      integer*4 now(3)
C
C --- Exibe horario do inicio da execuçao -----------------------------
C
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
      write (6, 400 )  (now(i),i=1,3)
  400 format ( 'Inicio ', i2, ':', i2, ':', i2)
C
C --- ABRE ARQUIVO ----------------------------------------------------
C
      OPEN (8, FILE='fig_laue.dat', STATUS='OLD', FORM='FORMATTED')
      OPEN (9, FILE='rede.dat', STATUS='NEW', FORM='FORMATTED')
C
C---- ENTRADA DE DADOS DAS MANCHAS DE LAUE ----------------------------
C
       M = 1
   10 READ(8,420,END=20) X(M),Y(M),ALIX
      write(6,*) M,X(M),Y(M)
  420 FORMAT(1X, 3F20.16)
      M = M+1
      GOTO 10
   20 CONTINUE
      CLOSE(8)
      MN = M - 1
C
C     NUMERO DE MANCHAS DE LAUE NO ANTEPARO: MN
C     NUMERO DE PONTOS NA GRADE: (2* NQC+1)**3 = NGRD ** 3 = NUM
C
      NQC = 50           ! Num ptos no quadrante da grade na posiçao do
C                          cristal
      ZL = 0.1           ! Dist cristal-anteparo em metros
      GRID = 0.2         ! Tamanho da grade no cristal em Angstrons
      ALBD = 0.67        ! Compr onda em Angstrons
C
      PI = 3.1415926538979D0
      D0 = 0.D+00
      D1 = 1.D+00
      D2 = 2.D+00
C
      NGRD = 2 * NQC + 1
      NUM = NGRD**3
      AK = D2 * PI / ALBD
      ZLQ = ZL * ZL
C
C --- Grava informaçoes da execuçao na primeira linha do arquivo ------
C
      WRITE(9,440) MN, ALBD,NUM
 440  FORMAT('#  Num manchas =',I7,4X,'Compr ond Angstr =',F6.3,4x,'Num
     *pts grade cristal = ',I7)
C
C---- CALCULO DA DENSIDADE RESULTANTE ---------------------------------
C
C ------ LOOPS PARA VARRER A GRADE NO CRISTAL -------------------------
C
        X2 = - DBLE(NQC + 1) * GRID
      DO 120 I = 1,NGRD
         X2 = X2 + GRID
C           X2 = -3.D+00
C      X2 = 0.D+00
      Y2 = - DBLE(NQC + 1) * GRID
      DO 130 J = 1,NGRD
         Y2 = Y2 + GRID
         Z2 = - DBLE(NQC + 1) * GRID
      DO 140 K = 1,NGRD
         Z2 = Z2 + GRID
C
      AMPC = D0
      AMPS = D0
C
C --- LOOP PARA VARRER AS MANCHAS DE LAUE NO ANTEPARO ----------------
C
      DO 160 L = 1,MN
C
        RR = DSQRT(X(L)**2 + Y(L)**2 + ZLQ)
        PROD = X(L)*X2 + Y(L)*Y2 + ZL*Z2
        FASE = AK * PROD / RR
        AMPC = AMPC + DCOS(FASE)
        AMPS = AMPS + DSIN(FASE)
C
  160 CONTINUE
C
        ACN = AMPC / DBLE(MN)     !Normaliza amplit result cosseno
        ACN2 = ACN * ACN
        ASN = AMPS / DBLE(MN)     !Normaliza ampl result seno
        ASN2 = ASN * ASN
        TOT = ACN2 + ASN2           !Intensidade total
C
C ----- Guarda em arq intensidades em cada ponto X2,Y2,Z2 do cristal ---
C
        WRITE(9,460) X2,Y2,Z2,ACN2,ASN2,TOT
  460   FORMAT (1X,6F20.16)
C
  140 CONTINUE
  130 CONTINUE
  120 CONTINUE
      CLOSE(9)
C
C --- Exibe horario do fim da execuçao e tempo gasto ------------------
C
      call itime(now)
      write (6, 480 )  (now(i), i=1,3)
  480 format ( 'Final ', i2, ':', i2, ':', i2)
      total = etime(elapsed)
      write(6,*) '#  tempo gasto', total
C
      END
