      SUBROUTINE CLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIST, TYPE;
      String             PATH;
      int                IMAT, KL, KU, M, MODE, N;
      REAL               ANORM, CNDNUM
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               SHRINK, TENTH
      const              SHRINK = 0.25E0, TENTH = 0.1E+0 ;
      REAL               ONE
      const              ONE = 1.0E+0 ;
      REAL               TWO
      const              TWO = 2.0E+0 ;
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

      // Set some constants for use in the subroutine.

      if ( FIRST ) {
         FIRST = .FALSE.
         EPS = SLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = SLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      }

      C2 = PATH( 2: 3 )

      // Set some parameters we don't plan to change.

      DIST = 'S'
      MODE = 3

      // xQR, xLQ, xQL, xRQ:  Set parameters to generate a general
                           // M x N matrix.

      if ( LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) {

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.1 ) {
            KL = 0
            KU = 0
         } else if ( IMAT.EQ.2 ) {
            KL = 0
            KU = MAX( N-1, 0 )
         } else if ( IMAT.EQ.3 ) {
            KL = MAX( M-1, 0 )
            KU = 0
         } else {
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         }

         // Set the condition number and norm.

         if ( IMAT.EQ.5 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.6 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.7 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.8 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'QK' ) ) {

         // xQK: truncated QR with pivoting.
              // Set parameters to generate a general
              // M x N matrix.

         // Set TYPE, the type of matrix to be generated.  'N' is nonsymmetric.

         TYPE = 'N'

         // Set DIST, the type of distribution for the random
         // number generator. 'S' is

         DIST = 'S'

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.2 ) {

            // 2. Random, Diagonal, CNDNUM = 2

            KL = 0
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         } else if ( IMAT.EQ.3 ) {

            // 3. Random, Upper triangular,  CNDNUM = 2

            KL = 0
            KU = MAX( N-1, 0 )
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         } else if ( IMAT.EQ.4 ) {

           // 4. Random, Lower triangular,  CNDNUM = 2

            KL = MAX( M-1, 0 )
            KU = 0
            CNDNUM = TWO
            ANORM = ONE
            MODE = 3
         } else {

            // 5.-19. Rectangular matrix

            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )

            if ( IMAT.GE.5 .AND. IMAT.LE.14 ) {

               // 5.-14. Random, CNDNUM = 2.

               CNDNUM = TWO
               ANORM = ONE
               MODE = 3

            } else if ( IMAT.EQ.15 ) {

               // 15. Random, CNDNUM = sqrt(0.1/EPS)

               CNDNUM = BADC1
               ANORM = ONE
               MODE = 3

            } else if ( IMAT.EQ.16 ) {

               // 16. Random, CNDNUM = 0.1/EPS

               CNDNUM = BADC2
               ANORM = ONE
               MODE = 3

            } else if ( IMAT.EQ.17 ) {

               // 17. Random, CNDNUM = 0.1/EPS,
                   // one small singular value S(N)=1/CNDNUM

               CNDNUM = BADC2
               ANORM = ONE
               MODE = 2

            } else if ( IMAT.EQ.18 ) {

               // 18. Random, scaled near underflow

               CNDNUM = TWO
               ANORM = SMALL
               MODE = 3

            } else if ( IMAT.EQ.19 ) {

               // 19. Random, scaled near overflow

               CNDNUM = TWO
               ANORM = LARGE
               MODE = 3

            }

         }

      } else if ( LSAMEN( 2, C2, 'GE' ) ) {

         // xGE:  Set parameters to generate a general M x N matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.1 ) {
            KL = 0
            KU = 0
         } else if ( IMAT.EQ.2 ) {
            KL = 0
            KU = MAX( N-1, 0 )
         } else if ( IMAT.EQ.3 ) {
            KL = MAX( M-1, 0 )
            KU = 0
         } else {
            KL = MAX( M-1, 0 )
            KU = MAX( N-1, 0 )
         }

         // Set the condition number and norm.

         if ( IMAT.EQ.8 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.9 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.10 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.11 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // xGB:  Set parameters to generate a general banded matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the condition number and norm.

         if ( IMAT.EQ.5 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.6 ) {
            CNDNUM = TENTH*BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.7 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.8 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'GT' ) ) {

         // xGT:  Set parameters to generate a general tridiagonal matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.1 ) {
            KL = 0
         } else {
            KL = 1
         }
         KU = KL

         // Set the condition number and norm.

         if ( IMAT.EQ.3 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.4 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.5 .OR. IMAT.EQ.11 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.6 .OR. IMAT.EQ.12 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'PP' ) ) {

         // xPO, xPP: Set parameters to generate a
         // symmetric or Hermitian positive definite matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = C2( 1: 1 )

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.1 ) {
            KL = 0
         } else {
            KL = MAX( N-1, 0 )
         }
         KU = KL

         // Set the condition number and norm.

         if ( IMAT.EQ.6 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.7 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.8 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.9 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'HE' ) .OR. LSAMEN( 2, C2, 'HP' ) .OR. LSAMEN( 2, C2, 'SY' ) .OR. LSAMEN( 2, C2, 'SP' ) ) {

         // xHE, xHP, xSY, xSP: Set parameters to generate a
         // symmetric or Hermitian matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = C2( 1: 1 )

         // Set the lower and upper bandwidths.

         if ( IMAT.EQ.1 ) {
            KL = 0
         } else {
            KL = MAX( N-1, 0 )
         }
         KU = KL

         // Set the condition number and norm.

         if ( IMAT.EQ.7 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.8 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.9 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.10 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // xPB:  Set parameters to generate a symmetric band matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'P'

         // Set the norm and condition number.

         if ( IMAT.EQ.5 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.6 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.7 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.8 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'PT' ) ) {

         // xPT:  Set parameters to generate a symmetric positive definite
        t // ridiagonal matrix.

         TYPE = 'P'
         if ( IMAT.EQ.1 ) {
            KL = 0
         } else {
            KL = 1
         }
         KU = KL

         // Set the condition number and norm.

         if ( IMAT.EQ.3 ) {
            CNDNUM = BADC1
         } else if ( IMAT.EQ.4 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( IMAT.EQ.5 .OR. IMAT.EQ.11 ) {
            ANORM = SMALL
         } else if ( IMAT.EQ.6 .OR. IMAT.EQ.12 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'TR' ) .OR. LSAMEN( 2, C2, 'TP' ) ) {

         // xTR, xTP:  Set parameters to generate a triangular matrix

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the lower and upper bandwidths.

         MAT = ABS( IMAT )
         if ( MAT.EQ.1 .OR. MAT.EQ.7 ) {
            KL = 0
            KU = 0
         } else if ( IMAT.LT.0 ) {
            KL = MAX( N-1, 0 )
            KU = 0
         } else {
            KL = 0
            KU = MAX( N-1, 0 )
         }

         // Set the condition number and norm.

         if ( MAT.EQ.3 .OR. MAT.EQ.9 ) {
            CNDNUM = BADC1
         } else if ( MAT.EQ.4 .OR. MAT.EQ.10 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( MAT.EQ.5 ) {
            ANORM = SMALL
         } else if ( MAT.EQ.6 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }

      } else if ( LSAMEN( 2, C2, 'TB' ) ) {

         // xTB:  Set parameters to generate a triangular band matrix.

         // Set TYPE, the type of matrix to be generated.

         TYPE = 'N'

         // Set the norm and condition number.

         MAT = ABS( IMAT )
         if ( MAT.EQ.2 .OR. MAT.EQ.8 ) {
            CNDNUM = BADC1
         } else if ( MAT.EQ.3 .OR. MAT.EQ.9 ) {
            CNDNUM = BADC2
         } else {
            CNDNUM = TWO
         }

         if ( MAT.EQ.4 ) {
            ANORM = SMALL
         } else if ( MAT.EQ.5 ) {
            ANORM = LARGE
         } else {
            ANORM = ONE
         }
      }
      IF( N.LE.1 ) CNDNUM = ONE

      RETURN

      // End of CLATB4

      }
