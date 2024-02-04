      void serrsy(PATH, NUNIT ) {

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
      REAL               ANRM, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX ), IW( NMAX );
      REAL               A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ), X( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, SSPCON, SSPRFS, SSPTRF, SSPTRI, SSPTRS, SSYCON, SSYCON_3, SSYCON_ROOK, SSYRFS, SSYTF2, SSYTF2_RK, SSYTF2_ROOK, SSYTRF, SSYTRF_RK, SSYTRF_ROOK, SSYTRI, SSYTRI_3, SSYTRI_3X, SSYTRI_ROOK, SSYTRI2, SSYTRI2X, SSYTRS, SSYTRS_3, SSYTRS_ROOK, SSYRFSX
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
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A[I, J] = 1. / REAL( I+J );
            AF[I, J] = 1. / REAL( I+J );
         } // 10
         B[J] = 0.0;
         E[J] = 0.0;
         R1[J] = 0.0;
         R2[J] = 0.0;
         W[J] = 0.0;
         X[J] = 0.0;
         IP[J] = J;
         IW[J] = J;
      } // 20
      ANRM = 1.0;
      RCOND = 1.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'SY' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with partial
         // (Bunch-Kaufman) pivoting.

         // SSYTRF

         SRNAMT = 'SSYTRF';
         INFOT = 1;
         ssytrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('SSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('SSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('SSYTRF', INFOT, NOUT, LERR, OK );

         // SSYTF2

         SRNAMT = 'SSYTF2';
         INFOT = 1;
         ssytf2('/', 0, A, 1, IP, INFO );
         chkxer('SSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytf2('U', -1, A, 1, IP, INFO );
         chkxer('SSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytf2('U', 2, A, 1, IP, INFO );
         chkxer('SSYTF2', INFOT, NOUT, LERR, OK );

         // SSYTRI

         SRNAMT = 'SSYTRI';
         INFOT = 1;
         ssytri('/', 0, A, 1, IP, W, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri('U', -1, A, 1, IP, W, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri('U', 2, A, 1, IP, W, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );

         // SSYTRI2

         SRNAMT = 'SSYTRI2';
         INFOT = 1;
         ssytri2('/', 0, A, 1, IP, W, IW, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri2('U', -1, A, 1, IP, W, IW, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri2('U', 2, A, 1, IP, W, IW, INFO );
         chkxer('SSYTRI', INFOT, NOUT, LERR, OK );

         // SSYTRI2X

         SRNAMT = 'SSYTRI2X';
         INFOT = 1;
         ssytri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRI2X', INFOT, NOUT, LERR, OK );

         // SSYTRS

         SRNAMT = 'SSYTRS';
         INFOT = 1;
         ssytrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('SSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('SSYTRS', INFOT, NOUT, LERR, OK );

         // SSYRFS

         SRNAMT = 'SSYRFS';
         INFOT = 1;
         ssyrfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyrfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssyrfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssyrfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssyrfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSYRFS', INFOT, NOUT, LERR, OK );

         // SSYRFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'SSYRFSX';
         INFOT = 1;
         ssyrfsx('/', EQ, 0, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssyrfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         EQ = 'N';
         INFOT = 3;
         ssyrfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssyrfsx('U', EQ, 0, -1, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssyrfsx('U', EQ, 2, 1, A, 1, AF, 2, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssyrfsx('U', EQ, 2, 1, A, 2, AF, 1, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         ssyrfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         ssyrfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('SSYRFSX', INFOT, NOUT, LERR, OK );

         // SSYCON

         SRNAMT = 'SSYCON';
         INFOT = 1;
         ssycon('/', 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssycon('U', -1, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssycon('U', 2, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssycon('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO );
         chkxer('SSYCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting.

         // SSYTRF_ROOK

         SRNAMT = 'SSYTRF_ROOK';
         INFOT = 1;
         ssytrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('SSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('SSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('SSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssytrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('SSYTRF_ROOK', INFOT, NOUT, LERR, OK );

         // SSYTF2_ROOK

         SRNAMT = 'SSYTF2_ROOK';
         INFOT = 1;
         ssytf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('SSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('SSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('SSYTF2_ROOK', INFOT, NOUT, LERR, OK );

         // SSYTRI_ROOK

         SRNAMT = 'SSYTRI_ROOK';
         INFOT = 1;
         ssytri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('SSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('SSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('SSYTRI_ROOK', INFOT, NOUT, LERR, OK );

         // SSYTRS_ROOK

         SRNAMT = 'SSYTRS_ROOK';
         INFOT = 1;
         ssytrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('SSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('SSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('SSYTRS_ROOK', INFOT, NOUT, LERR, OK );

         // SSYCON_ROOK

         SRNAMT = 'SSYCON_ROOK';
         INFOT = 1;
         ssycon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssycon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssycon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         ssycon_rook('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO );
         chkxer('SSYCON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // SSYTRF_RK

         SRNAMT = 'SSYTRF_RK';
         INFOT = 1;
         ssytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('SSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('SSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('SSYTRF_RK', INFOT, NOUT, LERR, OK );

         // SSYTF2_RK

         SRNAMT = 'SSYTF2_RK';
         INFOT = 1;
         ssytf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('SSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('SSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('SSYTF2_RK', INFOT, NOUT, LERR, OK );

         // SSYTRI_3

         SRNAMT = 'SSYTRI_3';
         INFOT = 1;
         ssytri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('SSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssytri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('SSYTRI_3', INFOT, NOUT, LERR, OK );

         // SSYTRI_3X

         SRNAMT = 'SSYTRI_3X';
         INFOT = 1;
         ssytri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssytri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('SSYTRI_3X', INFOT, NOUT, LERR, OK );

         // SSYTRS_3

         SRNAMT = 'SSYTRS_3';
         INFOT = 1;
         ssytrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('SSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssytrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('SSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssytrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('SSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         ssytrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('SSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         ssytrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('SSYTRS_3', INFOT, NOUT, LERR, OK );

         // SSYCON_3

         SRNAMT = 'SSYCON_3';
         INFOT = 1;
         ssycon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssycon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         ssycon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, IW, INFO);
         chkxer('SSYCON_3', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite packed matrix with partial
         // (Bunch-Kaufman) pivoting.

         // SSPTRF

         SRNAMT = 'SSPTRF';
         INFOT = 1;
         ssptrf('/', 0, A, IP, INFO );
         chkxer('SSPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssptrf('U', -1, A, IP, INFO );
         chkxer('SSPTRF', INFOT, NOUT, LERR, OK );

         // SSPTRI

         SRNAMT = 'SSPTRI';
         INFOT = 1;
         ssptri('/', 0, A, IP, W, INFO );
         chkxer('SSPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssptri('U', -1, A, IP, W, INFO );
         chkxer('SSPTRI', INFOT, NOUT, LERR, OK );

         // SSPTRS

         SRNAMT = 'SSPTRS';
         INFOT = 1;
         ssptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('SSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('SSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('SSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         ssptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('SSPTRS', INFOT, NOUT, LERR, OK );

         // SSPRFS

         SRNAMT = 'SSPRFS';
         INFOT = 1;
         ssprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         ssprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         ssprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         ssprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('SSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         ssprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('SSPRFS', INFOT, NOUT, LERR, OK );

         // SSPCON

         SRNAMT = 'SSPCON';
         INFOT = 1;
         sspcon('/', 0, A, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         sspcon('U', -1, A, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('SSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         sspcon('U', 1, A, IP, -1.0, RCOND, W, IW, INFO );
         chkxer('SSPCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
