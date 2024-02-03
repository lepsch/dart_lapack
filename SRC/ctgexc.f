      SUBROUTINE CTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, IFST, ILST, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               WANTQ, WANTZ;
      int                IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                HERE;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTGEX2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Decode and test input arguments.
      INFO = 0
      if ( N < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -7
      } else if ( LDQ < 1 || WANTQ && ( LDQ < MAX( 1, N ) ) ) {
         INFO = -9
      } else if ( LDZ < 1 || WANTZ && ( LDZ < MAX( 1, N ) ) ) {
         INFO = -11
      } else if ( IFST < 1 || IFST.GT.N ) {
         INFO = -12
      } else if ( ILST < 1 || ILST.GT.N ) {
         INFO = -13
      }
      if ( INFO != 0 ) {
         xerbla('CTGEXC', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.LE.1) RETURN       IF( IFST == ILST ) RETURN;

      if ( IFST < ILST ) {

         HERE = IFST

         } // 10

         // Swap with next one below

         ctgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO );
         if ( INFO != 0 ) {
            ILST = HERE
            RETURN
         }
         HERE = HERE + 1
         if (HERE < ILST) GO TO 10;
         HERE = HERE - 1
      } else {
         HERE = IFST - 1

         } // 20

         // Swap with next one above

         ctgex2(WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ, HERE, INFO );
         if ( INFO != 0 ) {
            ILST = HERE
            RETURN
         }
         HERE = HERE - 1
         if (HERE.GE.ILST) GO TO 20;
         HERE = HERE + 1
      }
      ILST = HERE
      RETURN

      // End of CTGEXC

      }
