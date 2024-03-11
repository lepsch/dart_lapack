      void zerrsy(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             PATH;
      int                infoc.NUNIT;
      // ..

      int                NMAX;
      const              NMAX = 4 ;
      String             C2;
      int                I, INFO, J;
      double             ANRM, RCOND;
      int                IP( NMAX );
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZSPCON, ZSPRFS, ZSPTRF, ZSPTRI, ZSPTRS, ZSYCON, ZSYCON_3, ZSYCON_ROOK, ZSYRFS, ZSYTF2, ZSYTF2_RK, ZSYTF2_ROOK, ZSYTRF, ZSYTRF_RK, ZSYTRF_ROOK, ZSYTRI, ZSYTRI_3, ZSYTRI_3X, ZSYTRI_ROOK, ZSYTRI2, ZSYTRI2X, ZSYTRS, ZSYTRS_3, ZSYTRS_ROOK
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, NOUT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX

      NOUT = infoc.NUNIT;
      NOUT.println( * );
      C2 = PATH.substring( 1, 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I][J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
         } // 10
         B[J] = 0.0;
         E[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         IP[J] = J;
      } // 20
      ANRM = 1.0;
      infoc.OK.value = true;

      if ( lsamen( 2, C2, 'SY' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // ZSYTRF

        srnamc.SRNAMT = 'ZSYTRF';
         infoc.INFOT = 1;
         zsytrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTF2

        srnamc.SRNAMT = 'ZSYTF2';
         infoc.INFOT = 1;
         zsytf2('/', 0, A, 1, IP, INFO );
         chkxer('ZSYTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytf2('U', -1, A, 1, IP, INFO );
         chkxer('ZSYTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytf2('U', 2, A, 1, IP, INFO );
         chkxer('ZSYTF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI

        srnamc.SRNAMT = 'ZSYTRI';
         infoc.INFOT = 1;
         zsytri('/', 0, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri('U', -1, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri('U', 2, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI2

        srnamc.SRNAMT = 'ZSYTRI2';
         infoc.INFOT = 1;
         zsytri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI2X

        srnamc.SRNAMT = 'ZSYTRI2X';
         infoc.INFOT = 1;
         zsytri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRS

        srnamc.SRNAMT = 'ZSYTRS';
         infoc.INFOT = 1;
         zsytrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsytrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsytrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZSYTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZSYTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYRFS

        srnamc.SRNAMT = 'ZSYRFS';
         infoc.INFOT = 1;
         zsyrfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsyrfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsyrfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsyrfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsyrfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYCON

        srnamc.SRNAMT = 'ZSYCON';
         infoc.INFOT = 1;
         zsycon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsycon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsycon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zsycon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SR' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // ZSYTRF_ROOK

        srnamc.SRNAMT = 'ZSYTRF_ROOK';
         infoc.INFOT = 1;
         zsytrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTF2_ROOK

        srnamc.SRNAMT = 'ZSYTF2_ROOK';
         infoc.INFOT = 1;
         zsytf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI_ROOK

        srnamc.SRNAMT = 'ZSYTRI_ROOK';
         infoc.INFOT = 1;
         zsytri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRS_ROOK

        srnamc.SRNAMT = 'ZSYTRS_ROOK';
         infoc.INFOT = 1;
         zsytrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsytrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsytrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZSYTRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYCON_ROOK

        srnamc.SRNAMT = 'ZSYCON_ROOK';
         infoc.INFOT = 1;
         zsycon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsycon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsycon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zsycon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // ZSYTRF_RK

        srnamc.SRNAMT = 'ZSYTRF_RK';
         infoc.INFOT = 1;
         zsytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('ZSYTRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZSYTRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZSYTRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTF2_RK

        srnamc.SRNAMT = 'ZSYTF2_RK';
         infoc.INFOT = 1;
         zsytf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI_3

        srnamc.SRNAMT = 'ZSYTRI_3';
         infoc.INFOT = 1;
         zsytri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZSYTRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZSYTRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRI_3X

        srnamc.SRNAMT = 'ZSYTRI_3X';
         infoc.INFOT = 1;
         zsytri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRS_3

        srnamc.SRNAMT = 'ZSYTRS_3';
         infoc.INFOT = 1;
         zsytrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsytrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsytrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('ZSYTRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zsytrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYCON_3

        srnamc.SRNAMT = 'ZSYCON_3';
         infoc.INFOT = 1;
         zsycon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsycon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsycon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, INFO);
         chkxer('ZSYCON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SP' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite packed matrix with partial
         // (Bunch-Kaufman) pivoting.

         // ZSPTRF

        srnamc.SRNAMT = 'ZSPTRF';
         infoc.INFOT = 1;
         zsptrf('/', 0, A, IP, INFO );
         chkxer('ZSPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsptrf('U', -1, A, IP, INFO );
         chkxer('ZSPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSPTRI

        srnamc.SRNAMT = 'ZSPTRI';
         infoc.INFOT = 1;
         zsptri('/', 0, A, IP, W, INFO );
         chkxer('ZSPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsptri('U', -1, A, IP, W, INFO );
         chkxer('ZSPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSPTRS

        srnamc.SRNAMT = 'ZSPTRS';
         infoc.INFOT = 1;
         zsptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSPRFS

        srnamc.SRNAMT = 'ZSPRFS';
         infoc.INFOT = 1;
         zsprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zsprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSPCON

        srnamc.SRNAMT = 'ZSPCON';
         infoc.INFOT = 1;
         zspcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zspcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zspcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'SA' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZSYTRF_AA

        srnamc.SRNAMT = 'ZSYTRF_AA';
         infoc.INFOT = 1;
         zsytrf_aa('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrf_aa('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytrf_aa('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf_aa('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrf_aa('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZSYTRS_AA

        srnamc.SRNAMT = 'ZSYTRS_AA';
         infoc.INFOT = 1;
         zsytrs_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrs_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsytrs_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsytrs_aa('U', 2, 1, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('ZSYTRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zsytrs_aa('U', 2, 1, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'S2' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZSYTRF_AA_2STAGE

        srnamc.SRNAMT = 'ZSYTRF_AA_2STAGE';
         infoc.INFOT = 1;
         zsytrf_aa_2stage('/', 0, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrf_aa_2stage('U', -1, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zsytrf_aa_2stage('U', 2, A, 1, A, 2, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zsytrf_aa_2stage('U', 2, A, 2, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zsytrf_aa_2stage('U', 2, A, 2, A, 8, IP, IP, W, 0, INFO );
         chkxer('ZSYTRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // CHETRS_AA_2STAGE

        srnamc.SRNAMT = 'ZSYTRS_AA_2STAGE';
         infoc.INFOT = 1;
         zsytrs_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zsytrs_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zsytrs_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zsytrs_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zsytrs_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zsytrs_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      }
