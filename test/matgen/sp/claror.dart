      void claror(SIDE, INIT, M, N, A, LDA, ISEED, X, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             INIT, SIDE;
      int                INFO, LDA, M, N;
      int                ISEED( 4 );
      Complex            A( LDA, * ), X( * );
      // ..

      double               ZERO, ONE, TOOSML;
      const              ZERO = 0.0, ONE = 1.0, TOOSML = 1.0e-20 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                IROW, ITYPE, IXFRM, J, JCOL, KBEG, NXFRM;
      double               FACTOR, XABS, XNORM;
      Complex            CSIGN, XNORMS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SCNRM2;
      //- COMPLEX            CLARND;
      // EXTERNAL lsame, SCNRM2, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CGERC, CLACGV, CLASET, CSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CMPLX, CONJG

      INFO = 0;
      if (N == 0 || M == 0) return;

      ITYPE = 0;
      if ( lsame( SIDE, 'L' ) ) {
         ITYPE = 1;
      } else if ( lsame( SIDE, 'R' ) ) {
         ITYPE = 2;
      } else if ( lsame( SIDE, 'C' ) ) {
         ITYPE = 3;
      } else if ( lsame( SIDE, 'T' ) ) {
         ITYPE = 4;
      }

      // Check for argument errors.

      if ( ITYPE == 0 ) {
         INFO = -1;
      } else if ( M < 0 ) {
         INFO = -3;
      } else if ( N < 0 || ( ITYPE == 3 && N != M ) ) {
         INFO = -4;
      } else if ( LDA < M ) {
         INFO = -6;
      }
      if ( INFO != 0 ) {
         xerbla('CLAROR', -INFO );
         return;
      }

      if ( ITYPE == 1 ) {
         NXFRM = M;
      } else {
         NXFRM = N;
      }

      // Initialize A to the identity matrix if desired

      if( lsame( INIT, 'I' ) ) claset( 'Full', M, N, CZERO, CONE, A, LDA );

      // If no rotation possible, still multiply by
      // a random complex number from the circle |x| = 1

       // 2)      Compute Rotation by computing Householder
               // Transformations H(2), H(3), ..., H(n).  Note that the
               // order in which they are computed is irrelevant.

      for (J = 1; J <= NXFRM; J++) { // 40
         X[J] = CZERO;
      } // 40

      for (IXFRM = 2; IXFRM <= NXFRM; IXFRM++) { // 60
         KBEG = NXFRM - IXFRM + 1;

         // Generate independent normal( 0, 1 ) random numbers

         for (J = KBEG; J <= NXFRM; J++) { // 50
            X[J] = CLARND( 3, ISEED );
         } // 50

         // Generate a Householder transformation from the random vector X

         XNORM = SCNRM2( IXFRM, X( KBEG ), 1 );
         XABS = ( X( KBEG ) ).abs();
         if ( XABS != CZERO ) {
            CSIGN = X( KBEG ) / XABS;
         } else {
            CSIGN = CONE;
         }
         XNORMS = CSIGN*XNORM;
         X[NXFRM+KBEG] = -CSIGN;
         FACTOR = XNORM*( XNORM+XABS );
         if ( ( FACTOR ).abs() < TOOSML ) {
            INFO = 1;
            xerbla('CLAROR', -INFO );
            return;
         } else {
            FACTOR = ONE / FACTOR;
         }
         X[KBEG] = X( KBEG ) + XNORMS;

         // Apply Householder transformation to A

         if ( ITYPE == 1 || ITYPE == 3 || ITYPE == 4 ) {

            // Apply H(k) on the left of A

            cgemv('C', IXFRM, N, CONE, A( KBEG, 1 ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            cgerc(IXFRM, N, -CMPLX( FACTOR ), X( KBEG ), 1, X( 2*NXFRM+1 ), 1, A( KBEG, 1 ), LDA );

         }

         if ( ITYPE >= 2 && ITYPE <= 4 ) {

            // Apply H(k)* (or H(k)') on the right of A

            if ( ITYPE == 4 ) {
               clacgv(IXFRM, X( KBEG ), 1 );
            }

            cgemv('N', M, IXFRM, CONE, A( 1, KBEG ), LDA, X( KBEG ), 1, CZERO, X( 2*NXFRM+1 ), 1 );
            cgerc(M, IXFRM, -CMPLX( FACTOR ), X( 2*NXFRM+1 ), 1, X( KBEG ), 1, A( 1, KBEG ), LDA );

         }
      } // 60

      X[1] = CLARND( 3, ISEED );
      XABS = ( X( 1 ) ).abs();
      if ( XABS != ZERO ) {
         CSIGN = X( 1 ) / XABS;
      } else {
         CSIGN = CONE;
      }
      X[2*NXFRM] = CSIGN;

      // Scale the matrix A by D.

      if ( ITYPE == 1 || ITYPE == 3 || ITYPE == 4 ) {
         for (IROW = 1; IROW <= M; IROW++) { // 70
            cscal(N, CONJG( X( NXFRM+IROW ) ), A( IROW, 1 ), LDA );
         } // 70
      }

      if ( ITYPE == 2 || ITYPE == 3 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 80
            cscal(M, X( NXFRM+JCOL ), A( 1, JCOL ), 1 );
         } // 80
      }

      if ( ITYPE == 4 ) {
         for (JCOL = 1; JCOL <= N; JCOL++) { // 90
            cscal(M, CONJG( X( NXFRM+JCOL ) ), A( 1, JCOL ), 1 );
         } // 90
      }
      return;
      }
