      void cerrhe(PATH, NUNIT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

// =====================================================================


      // .. Parameters ..
      int                NMAX;
      const              NMAX = 4 ;
      // ..
      // .. Local Scalars ..
      String             EQ;
      String             C2;
      int                I, INFO, J, N_ERR_BNDS, NPARAMS;
      double               ANRM, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double               R( NMAX ), R1( NMAX ), R2( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      Complex            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHECON, CHECON_3, CHECON_ROOK, CHERFS, CHETF2, CHETF2_RK, CHETF2_ROOK, CHETRF, CHETRF_RK, CHETRF_ROOK, CHETRI, CHETRI_3, CHETRI_3X, CHETRI_ROOK, CHETRI2, CHETRI2X, CHETRS, CHETRS_3, CHETRS_ROOK, CHKXER, CHPCON, CHPRFS, CHPTRF, CHPTRI, CHPTRS, CHERFSX
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NOUT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF[I, J] = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
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
      OK = true;

      if ( LSAMEN( 2, C2, 'HE' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // CHETRF

         SRNAMT = 'CHETRF';
         INFOT = 1;
         chetrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('CHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('CHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('CHETRF', INFOT, NOUT, LERR, OK );

         // CHETF2

         SRNAMT = 'CHETF2';
         INFOT = 1;
         chetf2('/', 0, A, 1, IP, INFO );
         chkxer('CHETF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetf2('U', -1, A, 1, IP, INFO );
         chkxer('CHETF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetf2('U', 2, A, 1, IP, INFO );
         chkxer('CHETF2', INFOT, NOUT, LERR, OK );

         // CHETRI

         SRNAMT = 'CHETRI';
         INFOT = 1;
         chetri('/', 0, A, 1, IP, W, INFO );
         chkxer('CHETRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri('U', -1, A, 1, IP, W, INFO );
         chkxer('CHETRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri('U', 2, A, 1, IP, W, INFO );
         chkxer('CHETRI', INFOT, NOUT, LERR, OK );

         // CHETRI2

         SRNAMT = 'CHETRI2';
         INFOT = 1;
         chetri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2', INFOT, NOUT, LERR, OK );

         // CHETRI2X

         SRNAMT = 'CHETRI2X';
         INFOT = 1;
         chetri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('CHETRI2X', INFOT, NOUT, LERR, OK );

         // CHETRS

         SRNAMT = 'CHETRS';
         INFOT = 1;
         chetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CHETRS', INFOT, NOUT, LERR, OK );

         // CHERFS

         SRNAMT = 'CHERFS';
         INFOT = 1;
         cherfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cherfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         cherfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cherfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         cherfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         cherfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cherfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CHERFS', INFOT, NOUT, LERR, OK );

         // CHECON

         SRNAMT = 'CHECON';
         INFOT = 1;
         checon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         checon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         checon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         checon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('CHECON', INFOT, NOUT, LERR, OK );

         // CHERFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'CHERFSX';
         INFOT = 1;
         cherfsx('/', EQ, 0, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cherfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         EQ = 'N';
         INFOT = 3;
         cherfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         cherfsx('U', EQ, 0, -1, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         cherfsx('U', EQ, 2, 1, A, 1, AF, 2, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         cherfsx('U', EQ, 2, 1, A, 2, AF, 1, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         cherfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         cherfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, R, INFO );
         chkxer('CHERFSX', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // CHETRF_ROOK

         SRNAMT = 'CHETRF_ROOK';
         INFOT = 1;
         chetrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('CHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('CHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chetrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('CHETRF_ROOK', INFOT, NOUT, LERR, OK );

         // CHETF2_ROOK

         SRNAMT = 'CHETF2_ROOK';
         INFOT = 1;
         chetf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('CHETF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('CHETF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('CHETF2_ROOK', INFOT, NOUT, LERR, OK );

         // CHETRI_ROOK

         SRNAMT = 'CHETRI_ROOK';
         INFOT = 1;
         chetri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('CHETRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('CHETRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('CHETRI_ROOK', INFOT, NOUT, LERR, OK );

         // CHETRS_ROOK

         SRNAMT = 'CHETRS_ROOK';
         INFOT = 1;
         chetrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CHETRS_ROOK', INFOT, NOUT, LERR, OK );

         // CHECON_ROOK

         SRNAMT = 'CHECON_ROOK';
         INFOT = 1;
         checon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         checon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         checon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         checon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('CHECON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HK' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // CHETRF_RK

         SRNAMT = 'CHETRF_RK';
         INFOT = 1;
         chetrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('CHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('CHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('CHETRF_RK', INFOT, NOUT, LERR, OK );

         // CHETF2_RK

         SRNAMT = 'CHETF2_RK';
         INFOT = 1;
         chetf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('CHETF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('CHETF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('CHETF2_RK', INFOT, NOUT, LERR, OK );

         // CHETRI_3

         SRNAMT = 'CHETRI_3';
         INFOT = 1;
         chetri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('CHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chetri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('CHETRI_3', INFOT, NOUT, LERR, OK );

         // CHETRI_3X

         SRNAMT = 'CHETRI_3X';
         INFOT = 1;
         chetri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         chetri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('CHETRI_3X', INFOT, NOUT, LERR, OK );

         // CHETRS_3

         SRNAMT = 'CHETRS_3';
         INFOT = 1;
         chetrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('CHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chetrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('CHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chetrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('CHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chetrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('CHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         chetrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('CHETRS_3', INFOT, NOUT, LERR, OK );

         // CHECON_3

         SRNAMT = 'CHECON_3';
         INFOT = 1;
         checon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         checon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         checon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('CHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         checon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, INFO);
         chkxer('CHECON_3', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

      // Test error exits of the routines that use factorization
      // of a Hermitian indefinite packed matrix with partial
      // (Bunch-Kaufman) diagonal pivoting method.

         // CHPTRF

         SRNAMT = 'CHPTRF';
         INFOT = 1;
         chptrf('/', 0, A, IP, INFO );
         chkxer('CHPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chptrf('U', -1, A, IP, INFO );
         chkxer('CHPTRF', INFOT, NOUT, LERR, OK );

         // CHPTRI

         SRNAMT = 'CHPTRI';
         INFOT = 1;
         chptri('/', 0, A, IP, W, INFO );
         chkxer('CHPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chptri('U', -1, A, IP, W, INFO );
         chkxer('CHPTRI', INFOT, NOUT, LERR, OK );

         // CHPTRS

         SRNAMT = 'CHPTRS';
         INFOT = 1;
         chptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('CHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('CHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('CHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         chptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('CHPTRS', INFOT, NOUT, LERR, OK );

         // CHPRFS

         SRNAMT = 'CHPRFS';
         INFOT = 1;
         chprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         chprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         chprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         chprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CHPRFS', INFOT, NOUT, LERR, OK );

         // CHPCON

         SRNAMT = 'CHPCON';
         INFOT = 1;
         chpcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('CHPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         chpcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('CHPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         chpcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('CHPCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
