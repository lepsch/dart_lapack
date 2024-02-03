      void derrsy(PATH, NUNIT ) {

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
      double             ANRM, RCOND, BERR;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX ), IW( NMAX );
      double             A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), R1( NMAX ), R2( NMAX ), W( 3*NMAX ), X( NMAX ), S( NMAX ), ERR_BNDS_N( NMAX, 3 ), ERR_BNDS_C( NMAX, 3 ), PARAMS( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, DSPCON, DSPRFS, DSPTRF, DSPTRI, DSPTRS, DSYCON, DSYCON_3, DSYCON_ROOK, DSYRFS, DSYTF2, DSYTF2_RK, DSYTF2_ROOK, DSYTRF, DSYTRF_RK, DSYTRF_ROOK, DSYTRI, DSYTRI_3, DSYTRI_3X, DSYTRI_ROOK, DSYTRI2, DSYTRI2X, DSYTRS, DSYTRS_3, DSYTRS_ROOK, DSYRFSX
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
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = 1.0 / DBLE( I+J );
            AF( I, J ) = 1.0 / DBLE( I+J );
         } // 10
         B( J ) = 0.0;
         E( J ) = 0.0;
         R1( J ) = 0.0;
         R2( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
         S( J ) = 0.0;
         IP( J ) = J;
         IW( J ) = J;
      } // 20
      ANRM = 1.0;
      RCOND = 1.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'SY' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with partial
         // (Bunch-Kaufman) pivoting.

         // DSYTRF

         SRNAMT = 'DSYTRF';
         INFOT = 1;
         dsytrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('DSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('DSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('DSYTRF', INFOT, NOUT, LERR, OK );

         // DSYTF2

         SRNAMT = 'DSYTF2';
         INFOT = 1;
         dsytf2('/', 0, A, 1, IP, INFO );
         chkxer('DSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytf2('U', -1, A, 1, IP, INFO );
         chkxer('DSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytf2('U', 2, A, 1, IP, INFO );
         chkxer('DSYTF2', INFOT, NOUT, LERR, OK );

         // DSYTRI

         SRNAMT = 'DSYTRI';
         INFOT = 1;
         dsytri('/', 0, A, 1, IP, W, INFO );
         chkxer('DSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri('U', -1, A, 1, IP, W, INFO );
         chkxer('DSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri('U', 2, A, 1, IP, W, INFO );
         chkxer('DSYTRI', INFOT, NOUT, LERR, OK );

         // DSYTRI2

         SRNAMT = 'DSYTRI2';
         INFOT = 1;
         dsytri2('/', 0, A, 1, IP, W, IW, INFO );
         chkxer('DSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri2('U', -1, A, 1, IP, W, IW, INFO );
         chkxer('DSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri2('U', 2, A, 1, IP, W, IW, INFO );
         chkxer('DSYTRI2', INFOT, NOUT, LERR, OK );

         // DSYTRI2X

         SRNAMT = 'DSYTRI2X';
         INFOT = 1;
         dsytri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRI2X', INFOT, NOUT, LERR, OK );

         // DSYTRS

         SRNAMT = 'DSYTRS';
         INFOT = 1;
         dsytrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('DSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('DSYTRS', INFOT, NOUT, LERR, OK );

         // DSYRFS

         SRNAMT = 'DSYRFS';
         INFOT = 1;
         dsyrfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyrfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsyrfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsyrfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsyrfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSYRFS', INFOT, NOUT, LERR, OK );

         // DSYRFSX

         N_ERR_BNDS = 3;
         NPARAMS = 0;
         SRNAMT = 'DSYRFSX';
         INFOT = 1;
         dsyrfsx('/', EQ, 0, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsyrfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         EQ = 'N';
         INFOT = 3;
         dsyrfsx('U', EQ, -1, 0, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsyrfsx('U', EQ, 0, -1, A, 1, AF, 1, IP, S, B, 1, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsyrfsx('U', EQ, 2, 1, A, 1, AF, 2, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsyrfsx('U', EQ, 2, 1, A, 2, AF, 1, IP, S, B, 2, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         dsyrfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 1, X, 2, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );
         INFOT = 14;
         dsyrfsx('U', EQ, 2, 1, A, 2, AF, 2, IP, S, B, 2, X, 1, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_N, ERR_BNDS_C, NPARAMS, PARAMS, W, IW, INFO );
         chkxer('DSYRFSX', INFOT, NOUT, LERR, OK );

         // DSYCON

         SRNAMT = 'DSYCON';
         INFOT = 1;
         dsycon('/', 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsycon('U', -1, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsycon('U', 2, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsycon('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO );
         chkxer('DSYCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting.

         // DSYTRF_ROOK

         SRNAMT = 'DSYTRF_ROOK';
         INFOT = 1;
         dsytrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('DSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('DSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('DSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsytrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('DSYTRF_ROOK', INFOT, NOUT, LERR, OK );

         // DSYTF2_ROOK

         SRNAMT = 'DSYTF2_ROOK';
         INFOT = 1;
         dsytf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('DSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('DSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('DSYTF2_ROOK', INFOT, NOUT, LERR, OK );

         // DSYTRI_ROOK

         SRNAMT = 'DSYTRI_ROOK';
         INFOT = 1;
         dsytri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('DSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('DSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('DSYTRI_ROOK', INFOT, NOUT, LERR, OK );

         // DSYTRS_ROOK

         SRNAMT = 'DSYTRS_ROOK';
         INFOT = 1;
         dsytrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('DSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('DSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('DSYTRS_ROOK', INFOT, NOUT, LERR, OK );

         // DSYCON_ROOK

         SRNAMT = 'DSYCON_ROOK';
         INFOT = 1;
         dsycon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsycon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsycon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         dsycon_rook('U', 1, A, 1, IP, -1.0, RCOND, W, IW, INFO);
         chkxer('DSYCON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // DSYTRF_RK

         SRNAMT = 'DSYTRF_RK';
         INFOT = 1;
         dsytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytrf_rk('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('DSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('DSYTRF_RK', INFOT, NOUT, LERR, OK );

         // DSYTF2_RK

         SRNAMT = 'DSYTF2_RK';
         INFOT = 1;
         dsytf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('DSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('DSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('DSYTF2_RK', INFOT, NOUT, LERR, OK );

         // DSYTRI_3

         SRNAMT = 'DSYTRI_3';
         INFOT = 1;
         dsytri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('DSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsytri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('DSYTRI_3', INFOT, NOUT, LERR, OK );

         // DSYTRI_3X

         SRNAMT = 'DSYTRI_3X';
         INFOT = 1;
         dsytri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsytri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('DSYTRI_3X', INFOT, NOUT, LERR, OK );

         // DSYTRS_3

         SRNAMT = 'DSYTRS_3';
         INFOT = 1;
         dsytrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('DSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsytrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('DSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsytrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('DSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dsytrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('DSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         dsytrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('DSYTRS_3', INFOT, NOUT, LERR, OK );

         // DSYCON_3

         SRNAMT = 'DSYCON_3';
         INFOT = 1;
         dsycon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsycon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         dsycon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, IW, INFO);
         chkxer('DSYCON_3', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite packed matrix with partial
         // (Bunch-Kaufman) pivoting.

         // DSPTRF

         SRNAMT = 'DSPTRF';
         INFOT = 1;
         dsptrf('/', 0, A, IP, INFO );
         chkxer('DSPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsptrf('U', -1, A, IP, INFO );
         chkxer('DSPTRF', INFOT, NOUT, LERR, OK );

         // DSPTRI

         SRNAMT = 'DSPTRI';
         INFOT = 1;
         dsptri('/', 0, A, IP, W, INFO );
         chkxer('DSPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsptri('U', -1, A, IP, W, INFO );
         chkxer('DSPTRI', INFOT, NOUT, LERR, OK );

         // DSPTRS

         SRNAMT = 'DSPTRS';
         INFOT = 1;
         dsptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('DSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('DSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('DSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         dsptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('DSPTRS', INFOT, NOUT, LERR, OK );

         // DSPRFS

         SRNAMT = 'DSPRFS';
         INFOT = 1;
         dsprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dsprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         dsprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         dsprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, IW, INFO );
         chkxer('DSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         dsprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, IW, INFO );
         chkxer('DSPRFS', INFOT, NOUT, LERR, OK );

         // DSPCON

         SRNAMT = 'DSPCON';
         INFOT = 1;
         dspcon('/', 0, A, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         dspcon('U', -1, A, IP, ANRM, RCOND, W, IW, INFO );
         chkxer('DSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         dspcon('U', 1, A, IP, -1.0, RCOND, W, IW, INFO );
         chkxer('DSPCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;
      }
