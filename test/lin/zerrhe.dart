      void zerrhe(PATH, infoc.NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                infoc.NUNIT;
      // ..

// =====================================================================


      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             C2;
      int                I, INFO, J;
      double             ANRM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      Complex         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZHECON, ZHECON_3, ZHECON_ROOK, ZHERFS, ZHETF2, ZHETF2_RK, ZHETF2_ROOK, ZHETRF, ZHETRF_RK, ZHETRF_ROOK, ZHETRF_AA, ZHETRF_AA_2STAGE, ZHETRI, ZHETRI_3, ZHETRI_3X, ZHETRI_ROOK, ZHETRI2, ZHETRI2X, ZHETRS, ZHETRS_3, ZHETRS_ROOK, ZHETRS_AA, ZHETRS_AA_2STAGE, ZHPCON, ZHPRFS, ZHPTRF, ZHPTRI, ZHPTRS
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
      // ..
      // .. Executable Statements ..

      NOUT = infoc.NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / (I+J).toDouble() );
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
      infoc.OK = true;

      if ( lsamen( 2, C2, 'HE' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // ZHETRF

        srnamc.SRNAMT = 'ZHETRF';
         infoc.INFOT = 1;
         zhetrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETF2

        srnamc.SRNAMT = 'ZHETF2';
         infoc.INFOT = 1;
         zhetf2('/', 0, A, 1, IP, INFO );
         chkxer('ZHETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetf2('U', -1, A, 1, IP, INFO );
         chkxer('ZHETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetf2('U', 2, A, 1, IP, INFO );
         chkxer('ZHETF2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI

        srnamc.SRNAMT = 'ZHETRI';
         infoc.INFOT = 1;
         zhetri('/', 0, A, 1, IP, W, INFO );
         chkxer('ZHETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri('U', -1, A, 1, IP, W, INFO );
         chkxer('ZHETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri('U', 2, A, 1, IP, W, INFO );
         chkxer('ZHETRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI2

        srnamc.SRNAMT = 'ZHETRI2';
         infoc.INFOT = 1;
         zhetri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI2X

        srnamc.SRNAMT = 'ZHETRI2X';
         infoc.INFOT = 1;
         zhetri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRS

        srnamc.SRNAMT = 'ZHETRS';
         infoc.INFOT = 1;
         zhetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhetrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhetrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZHETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZHETRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHERFS

        srnamc.SRNAMT = 'ZHERFS';
         infoc.INFOT = 1;
         zherfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zherfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zherfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zherfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zherfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zherfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 12;
         zherfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHECON

        srnamc.SRNAMT = 'ZHECON';
         infoc.INFOT = 1;
         zhecon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhecon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhecon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zhecon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHECON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'HR' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // ZHETRF_ROOK

        srnamc.SRNAMT = 'ZHETRF_ROOK';
         infoc.INFOT = 1;
         zhetrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETF2_ROOK

        srnamc.SRNAMT = 'ZHETF2_ROOK';
         infoc.INFOT = 1;
         zhetf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI_ROOK

        srnamc.SRNAMT = 'ZHETRI_ROOK';
         infoc.INFOT = 1;
         zhetri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRS_ROOK

        srnamc.SRNAMT = 'ZHETRS_ROOK';
         infoc.INFOT = 1;
         zhetrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhetrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhetrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZHETRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHECON_ROOK

        srnamc.SRNAMT = 'ZHECON_ROOK';
         infoc.INFOT = 1;
         zhecon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhecon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhecon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zhecon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'HK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // ZHETRF_RK

        srnamc.SRNAMT = 'ZHETRF_RK';
         infoc.INFOT = 1;
         zhetrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('ZHETRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZHETRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZHETRF_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETF2_RK

        srnamc.SRNAMT = 'ZHETF2_RK';
         infoc.INFOT = 1;
         zhetf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI_3

        srnamc.SRNAMT = 'ZHETRI_3';
         infoc.INFOT = 1;
         zhetri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZHETRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZHETRI_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRI_3X

        srnamc.SRNAMT = 'ZHETRI_3X';
         infoc.INFOT = 1;
         zhetri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRS_3

        srnamc.SRNAMT = 'ZHETRS_3';
         infoc.INFOT = 1;
         zhetrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhetrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhetrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('ZHETRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 9;
         zhetrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHECON_3

        srnamc.SRNAMT = 'ZHECON_3';
         infoc.INFOT = 1;
         zhecon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhecon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhecon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhecon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, INFO);
         chkxer('ZHECON_3', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with Aasen's algorithm.

      } else if ( lsamen( 2, C2, 'HA' ) ) {

         // ZHETRF_AA

        srnamc.SRNAMT = 'ZHETRF_AA';
         infoc.INFOT = 1;
         zhetrf_aa('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrf_aa('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetrf_aa('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf_aa('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrf_aa('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRS_AA

        srnamc.SRNAMT = 'ZHETRS_AA';
         infoc.INFOT = 1;
         zhetrs_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrs_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhetrs_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhetrs_aa('U', 2, 1, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhetrs_aa('U', 2, 1, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zhetrs_aa('U', 0, 1, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zhetrs_aa('U', 0, 1, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('ZHETRS_AA', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'S2' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZHETRF_AA_2STAGE

        srnamc.SRNAMT = 'ZHETRF_AA_2STAGE';
         infoc.INFOT = 1;
         zhetrf_aa_2stage('/', 0, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrf_aa_2stage('U', -1, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 4;
         zhetrf_aa_2stage('U', 2, A, 1, A, 2, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 6;
         zhetrf_aa_2stage('U', 2, A, 2, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zhetrf_aa_2stage('U', 2, A, 2, A, 8, IP, IP, W, 0, INFO );
         chkxer('ZHETRF_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHETRS_AA_2STAGE

        srnamc.SRNAMT = 'ZHETRS_AA_2STAGE';
         infoc.INFOT = 1;
         zhetrs_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhetrs_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhetrs_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhetrs_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhetrs_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 11;
         zhetrs_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_STAGE', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

      } else if ( lsamen( 2, C2, 'HP' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite packed matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // ZHPTRF

        srnamc.SRNAMT = 'ZHPTRF';
         infoc.INFOT = 1;
         zhptrf('/', 0, A, IP, INFO );
         chkxer('ZHPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhptrf('U', -1, A, IP, INFO );
         chkxer('ZHPTRF', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHPTRI

        srnamc.SRNAMT = 'ZHPTRI';
         infoc.INFOT = 1;
         zhptri('/', 0, A, IP, W, INFO );
         chkxer('ZHPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhptri('U', -1, A, IP, W, INFO );
         chkxer('ZHPTRI', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHPTRS

        srnamc.SRNAMT = 'ZHPTRS';
         infoc.INFOT = 1;
         zhptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 7;
         zhptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHPRFS

        srnamc.SRNAMT = 'ZHPRFS';
         infoc.INFOT = 1;
         zhprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 3;
         zhprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 8;
         zhprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 10;
         zhprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );

         // ZHPCON

        srnamc.SRNAMT = 'ZHPCON';
         infoc.INFOT = 1;
         zhpcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 2;
         zhpcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
         infoc.INFOT = 5;
         zhpcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', infoc.INFOT, NOUT, infoc.LERR, infoc.OK );
      }

      // Print a summary line.

      alaesm(PATH, infoc.OK, NOUT );

      return;
      }
