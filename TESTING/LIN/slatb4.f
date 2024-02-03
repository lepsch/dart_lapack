      SUBROUTINE SLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             DIST, TYPE;
      String             PATH;
      int                IMAT, KL, KU, M, MODE, N;
      REAL               ANORM, CNDNUM
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               SHRINK, TENTH
      PARAMETER          ( SHRINK = 0.25E0, TENTH = 0.1E+0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               FIRST;
      String             C2;
      int                MAT;
      REAL               BADC1, BADC2, EPS, LARGE, SMALL
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      REAL               SLAMCH
      // EXTERNAL LSAMEN, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
      // ..
      // .. Data statements ..
      DATA               FIRST / .TRUE. /
      // ..
      // .. Executable Statements ..
*
      // Set some constants for use in the subroutine.
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         EPS = SLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = SLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      END IF
*
      C2 = PATH( 2: 3 )
*
      // Set some parameters we don't plan to change.
*
      DIST = 'S'
      MODE = 3
*
      IF( LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) THEN
*
         // xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
                              // M x N matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.EQ.2 ) THEN
            KL = 0
            KU = MAX( N-1, 0 )
         ELSE IF( IMAT.EQ.3 ) THEN
            KL = MAX( M-1, 0 )
            KU = 0
         ELSE
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         END IF
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'QK' ) ) THEN
*
         // xQK: truncated QR with pivoting.
              // Set parameters to generate a general
              // M x N matrix.
*
         // Set TYPE, the type of matrix to be generated.  'N' is nonsymmetric.
*
         TYPE = 'N'
*
         // Set DIST, the type of distribution for the random
         // number generator. 'S' is
*
         DIST = 'S'
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.2 ) THEN
*
            // 2. Random, Diagonal, CNDNUM = 2
*
            KL = 0
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE IF( IMAT.EQ.3 ) THEN
*
            // 3. Random, Upper triangular,  CNDNUM = 2
*
            KL = 0
            KU = MAX( N-1, 0 )
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE IF( IMAT.EQ.4 ) THEN
*
           // 4. Random, Lower triangular,  CNDNUM = 2
*
            KL = MAX( M-1, 0 )
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         ELSE
*
            // 5.-19. Rectangular matrix
*
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
*
            IF( IMAT.GE.5 .AND. IMAT.LE.14 ) THEN
*
               // 5.-14. Random, CNDNUM = 2.
*
               CNDNUM = TWO
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.15 ) THEN
*
               // 15. Random, CNDNUM = sqrt(0.1/EPS)
*
               CNDNUM = BADC1
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.16 ) THEN
*
               // 16. Random, CNDNUM = 0.1/EPS
*
               CNDNUM = BADC2
               ANORM = ONE
               MODE = 3
*
            ELSE IF( IMAT.EQ.17 ) THEN
*
               // 17. Random, CNDNUM = 0.1/EPS,
                   // one small singular value S(N)=1/CNDNUM
*
               CNDNUM = BADC2
               ANORM = ONE
               MODE = 2
*
            ELSE IF( IMAT.EQ.18 ) THEN
*
               // 18. Random, scaled near underflow
*
               CNDNUM = TWO
               ANORM = SMALL
               MODE = 3
*
            ELSE IF( IMAT.EQ.19 ) THEN
*
               // 19. Random, scaled near overflow
*
               CNDNUM = TWO
               ANORM = LARGE
               MODE = 3
*
            END IF
*
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN
*
         // xGE:  Set parameters to generate a general M x N matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.EQ.2 ) THEN
            KL = 0
            KU = MAX( N-1, 0 )
         ELSE IF( IMAT.EQ.3 ) THEN
            KL = MAX( M-1, 0 )
            KU = 0
         ELSE
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         END IF
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.8 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.9 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.10 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.11 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
         // xGB:  Set parameters to generate a general banded matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = TENTH*BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN
*
         // xGT:  Set parameters to generate a general tridiagonal matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = 1
         END IF
         KU = KL
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.3 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.5 .OR. IMAT.EQ.11 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.6 .OR. IMAT.EQ.12 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'PP' ) ) THEN
*
         // xPO, xPP, xSY, xSP: Set parameters to generate a
         // symmetric positive definite matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = C2( 1: 1 )
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = MAX( N-1, 0 )
         END IF
         KU = KL
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.7 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.8 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.9 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
*
      ELSE IF( LSAMEN( 2, C2, 'SY' ) .OR. LSAMEN( 2, C2, 'SP' ) ) THEN
*
         // xSY, xSP: Set parameters to generate a
         // symmetric matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = C2( 1: 1 )
*
         // Set the lower and upper bandwidths.
*
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = MAX( N-1, 0 )
         END IF
         KU = KL
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.7 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.8 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.9 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.10 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
*
         // xPB:  Set parameters to generate a symmetric band matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'P'
*
         // Set the norm and condition number.
*
         IF( IMAT.EQ.5 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN
*
         // xPT:  Set parameters to generate a symmetric positive definite
        t // ridiagonal matrix.
*
         TYPE = 'P'
         IF( IMAT.EQ.1 ) THEN
            KL = 0
         ELSE
            KL = 1
         END IF
         KU = KL
*
         // Set the condition number and norm.
*
         IF( IMAT.EQ.3 ) THEN
            CNDNUM = BADC1
         ELSE IF( IMAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( IMAT.EQ.5 .OR. IMAT.EQ.11 ) THEN
            ANORM = SMALL
         ELSE IF( IMAT.EQ.6 .OR. IMAT.EQ.12 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) .OR. LSAMEN( 2, C2, 'TP' ) ) THEN
*
         // xTR, xTP:  Set parameters to generate a triangular matrix
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the lower and upper bandwidths.
*
         MAT = ABS( IMAT )
         IF( MAT.EQ.1 .OR. MAT.EQ.7 ) THEN
            KL = 0
            KU = 0
         ELSE IF( IMAT.LT.0 ) THEN
            KL = MAX( N-1, 0 )
            KU = 0
         ELSE
            KL = 0
            KU = MAX( N-1, 0 )
         END IF
*
         // Set the condition number and norm.
*
         IF( MAT.EQ.3 .OR. MAT.EQ.9 ) THEN
            CNDNUM = BADC1
         ELSE IF( MAT.EQ.4 ) THEN
            CNDNUM = BADC2
         ELSE IF( MAT.EQ.10 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( MAT.EQ.5 ) THEN
            ANORM = SMALL
         ELSE IF( MAT.EQ.6 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
         // xTB:  Set parameters to generate a triangular band matrix.
*
         // Set TYPE, the type of matrix to be generated.
*
         TYPE = 'N'
*
         // Set the norm and condition number.
*
         MAT = ABS( IMAT )
         IF( MAT.EQ.2 .OR. MAT.EQ.8 ) THEN
            CNDNUM = BADC1
         ELSE IF( MAT.EQ.3 .OR. MAT.EQ.9 ) THEN
            CNDNUM = BADC2
         ELSE
            CNDNUM = TWO
         END IF
*
         IF( MAT.EQ.4 ) THEN
            ANORM = SMALL
         ELSE IF( MAT.EQ.5 ) THEN
            ANORM = LARGE
         ELSE
            ANORM = ONE
         END IF
      END IF
      IF( N.LE.1 ) CNDNUM = ONE
*
      RETURN
*
      // End of SLATB4
*
      END
