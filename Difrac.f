      PROGRAM DIFRAC
C
C-----------------------------------------------------------------------
C     SIMULACAO COMPUTACIONAL DA DIFRACAO DE LAUE 
C     A PARTIR DA REDE CRISTALINA SC
C
C     Calcula intensidades numa grade do anteparo
C
C     9: OUTPUT Intensidade para cada X2,Y2 do anteparo 'out.dat'
C
C     Osmar S.Silva Jr. , 2017/08/29
C-----------------------------------------------------------------------
C
C     DIMENSION X(200)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ REDE/ X(200),NQ,NAT,DIST
C
C --- ABRE ARQUIVOS ---------------------------------------------------
C
      OPEN (8, FILE='in.dat', STATUS='OLD', FORM='FORMATTED')
      OPEN (9, FILE='out.dat', STATUS='NEW', FORM='FORMATTED')
C
C---- ENTRADA DE DADOS ------------------------------------------------
C
      READ(8,410) NQ,NQA,Z2,GRID
  410 FORMAT(5X,I4,5X,I7,4X,F20.16,6X,F20.16)
      READ(8,415) DIST,ALBD
  415 FORMAT(7X,F20.16,6X,F20.16)
C
C     NUMERO DE ATOMOS: (2* NQ +1)**3 = NAT**3 = NAT3
C     NUMERO DE PONTOS NO ANTEPARO: (2* NQA+1)**2 = NANT*NANT = NANT2
C
C     NQ = 20
C     NQA = 100
C     Z2 = 1.D-01        !Dist cristal-anteparo em metros
C     GRID = 1.D-03      !Tamanho da grade no anteparo em metros
C     DIST = 3.35D+00    !Dist atomo-atomo em Angstrons
C     ALBD = 0.67D+00    !Compr onda em Angstrons
C
      PI = 3.1415926538979D0
      D0 = 0.D+00
      D2 = 2.D+00
      NAT = 2 * NQ + 1
      NAT3 = NAT ** 3
      NANT = 2 * NQA + 1
      NANT2 = NANT * NANT
      AK = D2 * PI / ALBD
      Z2Q = Z2 * Z2
C
C --- Grava informaçoes da execuçao na primeira linha do arquivo ------
C
      WRITE(9,420) NQ,NQA,Z2,GRID
  420 FORMAT('# NQ=',I4,' NQA=',I7,' Z2=',F20.16,' GRID=',F20.16)
      WRITE(9,425) DIST,ALBD
  425 FORMAT('# DIST=',F20.16,' ALBD=',F20.16)
C
C---- DEFININDO COORDENADAS DOS ATOMOS --------------------------------
C
      CALL SC
C
C---- CALCULO DA AMPLITUDE RESULTANTE ---------------------------------
C
        X2 = - DBLE(NQA + 1) * GRID
C
C --- LOOPS PARA VARRER A GRADE NO ANTEPARO -------------------------
C
      DO 100 IX2 = 1,NANT
        X2 = X2 + GRID
        Y2 = - DBLE(NQA + 1) * GRID
C
      DO 110 IY2 = 1,NANT
        Y2 = Y2 + GRID
        AMPC = D0
        AMPS = D0
        RR = DSQRT(X2**2 + Y2**2 + Z2Q)
C
C --- LOOPS PARA COBRIR TODOS OS ATOMOS DO CRISTAL --------------------
C
      DO 120 I = 1,NAT
      DO 130 J = 1,NAT
      DO 140 K = 1,NAT
C
        PROD = X(I)*X2+X(J)*Y2+X(K)*Z2
        FASE = AK * PROD / RR
        AMPC = AMPC + DCOS(FASE)
        AMPS = AMPS + DSIN(FASE)
C
  140 CONTINUE
  130 CONTINUE
  120 CONTINUE
C
        ACN = AMPC / DBLE(NAT3)     !Normaliza amplit result cosseno
        ACN2 = ACN * ACN
        ASN = AMPS / DBLE(NAT3)     !Normaliza ampl result seno
        ASN2 = ASN * ASN
        TOT = ACN2 + ASN2           !Intensidade total
C
C ----- Guarda em arq intensidades em cada ponto X2,Y2 do anteparo ----
C
        WRITE(9,460) X2,Y2,ACN2,ASN2,TOT
  460   FORMAT (1X,5F20.16)
C
  110 CONTINUE
  100 CONTINUE
C
      CLOSE(9)
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SC
C
C     COORDENADAS DOS ATOMOS PARA A REDE CRISTALINA CUBICA SIMPLES
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ REDE/ X(200),NQ,NAT,DIST
C
      X(1) = - DBLE(NQ) * DIST
      DO 10 I = 2,NAT
        X(I) = X(I-1) + DIST
  10  CONTINUE
C
      RETURN
      END
