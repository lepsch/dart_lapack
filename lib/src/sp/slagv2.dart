      void slagv2(final int A, final int LDA, final Matrix<double> B, final int LDB, final int ALPHAR, final int ALPHAI, final int BETA, final int CSL, final int SNL, final int CSR, final int SNR,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB;
      double               CSL, CSR, SNL, SNR;
      double               A( LDA, * ), ALPHAI( 2 ), ALPHAR( 2 ), B( LDB, * ), BETA( 2 );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ANORM, ASCALE, BNORM, BSCALE, H1, H2, H3, QQ, R, RR, SAFMIN, SCALE1, SCALE2, T, ULP, WI, WR1, WR2;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAG2, SLARTG, SLASV2, SROT
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLAPY2;
      // EXTERNAL SLAMCH, SLAPY2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX

      SAFMIN = SLAMCH( 'S' );
      ULP = SLAMCH( 'P' );

      // Scale A

      ANORM = max( ( A( 1, 1 ) ).abs()+( A( 2, 1 ) ).abs(), ( A( 1, 2 ) ).abs()+( A( 2, 2 ) ).abs(), SAFMIN );
      ASCALE = ONE / ANORM;
      A[1][1] = ASCALE*A( 1, 1 );
      A[1][2] = ASCALE*A( 1, 2 );
      A[2][1] = ASCALE*A( 2, 1 );
      A[2][2] = ASCALE*A( 2, 2 );

      // Scale B

      BNORM = max( ( B( 1, 1 ) ).abs(), ( B( 1, 2 ) ).abs()+( B( 2, 2 ) ).abs(), SAFMIN );
      BSCALE = ONE / BNORM;
      B[1][1] = BSCALE*B( 1, 1 );
      B[1][2] = BSCALE*B( 1, 2 );
      B[2][2] = BSCALE*B( 2, 2 );

      // Check if A can be deflated

      if ( ( A( 2, 1 ) ).abs() <= ULP ) {
         CSL = ONE;
         SNL = ZERO;
         CSR = ONE;
         SNR = ZERO;
         A[2][1] = ZERO;
         B[2][1] = ZERO;
         WI = ZERO;

      // Check if B is singular

      } else if ( ( B( 1, 1 ) ).abs() <= ULP ) {
         slartg(A( 1, 1 ), A( 2, 1 ), CSL, SNL, R );
         CSR = ONE;
         SNR = ZERO;
         srot(2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL );
         srot(2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL );
         A[2][1] = ZERO;
         B[1][1] = ZERO;
         B[2][1] = ZERO;
         WI = ZERO;

      } else if ( ( B( 2, 2 ) ).abs() <= ULP ) {
         slartg(A( 2, 2 ), A( 2, 1 ), CSR, SNR, T );
         SNR = -SNR;
         srot(2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR );
         srot(2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR );
         CSL = ONE;
         SNL = ZERO;
         A[2][1] = ZERO;
         B[2][1] = ZERO;
         B[2][2] = ZERO;
         WI = ZERO;

      } else {

         // B is nonsingular, first compute the eigenvalues of (A,B)

         slag2(A, LDA, B, LDB, SAFMIN, SCALE1, SCALE2, WR1, WR2, WI );

         if ( WI == ZERO ) {

            // two real eigenvalues, compute s*A-w*B

            H1 = SCALE1*A( 1, 1 ) - WR1*B( 1, 1 );
            H2 = SCALE1*A( 1, 2 ) - WR1*B( 1, 2 );
            H3 = SCALE1*A( 2, 2 ) - WR1*B( 2, 2 );

            RR = SLAPY2( H1, H2 );
            QQ = SLAPY2( SCALE1*A( 2, 1 ), H3 );

            if ( RR > QQ ) {

               // find right rotation matrix to zero 1,1 element of
               // (sA - wB)

               slartg(H2, H1, CSR, SNR, T );

            } else {

               // find right rotation matrix to zero 2,1 element of
               // (sA - wB)

               slartg(H3, SCALE1*A( 2, 1 ), CSR, SNR, T );

            }

            SNR = -SNR;
            srot(2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR );
            srot(2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR );

            // compute inf norms of A and B

            H1 = max( ( A( 1, 1 ) ).abs()+( A( 1, 2 ) ).abs(), ( A( 2, 1 ) ).abs()+( A( 2, 2 ) ).abs() )             H2 = max( ( B( 1, 1 ) ).abs()+( B( 1, 2 ) ).abs(), ( B( 2, 1 ) ).abs()+( B( 2, 2 ) ).abs() );

            if ( ( SCALE1*H1 ) >= ( WR1 ).abs()*H2 ) {

               // find left rotation matrix Q to zero out B(2,1)

               slartg(B( 1, 1 ), B( 2, 1 ), CSL, SNL, R );

            } else {

               // find left rotation matrix Q to zero out A(2,1)

               slartg(A( 1, 1 ), A( 2, 1 ), CSL, SNL, R );

            }

            srot(2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL );
            srot(2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL );

            A[2][1] = ZERO;
            B[2][1] = ZERO;

         } else {

            // a pair of complex conjugate eigenvalues
            // first compute the SVD of the matrix B

            slasv2(B( 1, 1 ), B( 1, 2 ), B( 2, 2 ), R, T, SNR, CSR, SNL, CSL );

            // Form (A,B) := Q(A,B)Z**T where Q is left rotation matrix and
            // Z is right rotation matrix computed from SLASV2

            srot(2, A( 1, 1 ), LDA, A( 2, 1 ), LDA, CSL, SNL );
            srot(2, B( 1, 1 ), LDB, B( 2, 1 ), LDB, CSL, SNL );
            srot(2, A( 1, 1 ), 1, A( 1, 2 ), 1, CSR, SNR );
            srot(2, B( 1, 1 ), 1, B( 1, 2 ), 1, CSR, SNR );

            B[2][1] = ZERO;
            B[1][2] = ZERO;

         }

      }

      // Unscaling

      A[1][1] = ANORM*A( 1, 1 );
      A[2][1] = ANORM*A( 2, 1 );
      A[1][2] = ANORM*A( 1, 2 );
      A[2][2] = ANORM*A( 2, 2 );
      B[1][1] = BNORM*B( 1, 1 );
      B[2][1] = BNORM*B( 2, 1 );
      B[1][2] = BNORM*B( 1, 2 );
      B[2][2] = BNORM*B( 2, 2 );

      if ( WI == ZERO ) {
         ALPHAR[1] = A( 1, 1 );
         ALPHAR[2] = A( 2, 2 );
         ALPHAI[1] = ZERO;
         ALPHAI[2] = ZERO;
         BETA[1] = B( 1, 1 );
         BETA[2] = B( 2, 2 );
      } else {
         ALPHAR[1] = ANORM*WR1 / SCALE1 / BNORM;
         ALPHAI[1] = ANORM*WI / SCALE1 / BNORM;
         ALPHAR[2] = ALPHAR( 1 );
         ALPHAI[2] = -ALPHAI( 1 );
         BETA[1] = ONE;
         BETA[2] = ONE;
      }

      }
