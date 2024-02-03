      SUBROUTINE CERRSY( PATH, NUNIT );

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
      String             C2;
      int                I, INFO, J;
      REAL               ANRM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      REAL               R( NMAX ), R1( NMAX ), R2( NMAX );
      COMPLEX            A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, CSPCON, CSPRFS, CSPTRF, CSPTRI, CSPTRS, CSYCON, CSYCON_3, CSYCON_ROOK, CSYRFS, CSYTF2, CSYTF2_RK, CSYTF2_ROOK, CSYTRF, CSYTRF_RK, CSYTRF_ROOK, CSYTRI, CSYTRI_3, CSYTRI_3X, CSYTRI_ROOK, CSYTRI2, CSYTRI2X, CSYTRS, CSYTRS_3, CSYTRS_ROOK
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
            A( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
            AF( I, J ) = CMPLX( 1. / REAL( I+J ), -1. / REAL( I+J ) );
         } // 10
         B( J ) = 0.0;
         E( J ) = 0.0;
         R1( J ) = 0.0;
         R2( J ) = 0.0;
         W( J ) = 0.0;
         X( J ) = 0.0;
         IP( J ) = J;
      } // 20
      ANRM = 1.0;
      OK = true;

      if ( LSAMEN( 2, C2, 'SY' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // CSYTRF

         SRNAMT = 'CSYTRF';
         INFOT = 1;
         csytrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('CSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('CSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('CSYTRF', INFOT, NOUT, LERR, OK );

         // CSYTF2

         SRNAMT = 'CSYTF2';
         INFOT = 1;
         csytf2('/', 0, A, 1, IP, INFO );
         chkxer('CSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytf2('U', -1, A, 1, IP, INFO );
         chkxer('CSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytf2('U', 2, A, 1, IP, INFO );
         chkxer('CSYTF2', INFOT, NOUT, LERR, OK );

         // CSYTRI

         SRNAMT = 'CSYTRI';
         INFOT = 1;
         csytri('/', 0, A, 1, IP, W, INFO );
         chkxer('CSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri('U', -1, A, 1, IP, W, INFO );
         chkxer('CSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri('U', 2, A, 1, IP, W, INFO );
         chkxer('CSYTRI', INFOT, NOUT, LERR, OK );

         // CSYTRI2

         SRNAMT = 'CSYTRI2';
         INFOT = 1;
         csytri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2', INFOT, NOUT, LERR, OK );

         // CSYTRI2X

         SRNAMT = 'CSYTRI2X';
         INFOT = 1;
         csytri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRI2X', INFOT, NOUT, LERR, OK );

         // CSYTRS

         SRNAMT = 'CSYTRS';
         INFOT = 1;
         csytrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csytrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csytrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CSYTRS', INFOT, NOUT, LERR, OK );

         // CSYRFS

         SRNAMT = 'CSYRFS';
         INFOT = 1;
         csyrfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csyrfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csyrfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csyrfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csyrfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         csyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CSYRFS', INFOT, NOUT, LERR, OK );

         // CSYCON

         SRNAMT = 'CSYCON';
         INFOT = 1;
         csycon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csycon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csycon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         csycon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('CSYCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // CSYTRF_ROOK

         SRNAMT = 'CSYTRF_ROOK';
         INFOT = 1;
         csytrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('CSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('CSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('CSYTRF_ROOK', INFOT, NOUT, LERR, OK );

         // CSYTF2_ROOK

         SRNAMT = 'CSYTF2_ROOK';
         INFOT = 1;
         csytf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('CSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('CSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('CSYTF2_ROOK', INFOT, NOUT, LERR, OK );

         // CSYTRI_ROOK

         SRNAMT = 'CSYTRI_ROOK';
         INFOT = 1;
         csytri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('CSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('CSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('CSYTRI_ROOK', INFOT, NOUT, LERR, OK );

         // CSYTRS_ROOK

         SRNAMT = 'CSYTRS_ROOK';
         INFOT = 1;
         csytrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csytrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('CSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csytrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('CSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('CSYTRS_ROOK', INFOT, NOUT, LERR, OK );

         // CSYCON_ROOK

         SRNAMT = 'CSYCON_ROOK';
         INFOT = 1;
         csycon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csycon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csycon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         csycon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('CSYCON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // CSYTRF_RK

         SRNAMT = 'CSYTRF_RK';
         INFOT = 1;
         csytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('CSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('CSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('CSYTRF_RK', INFOT, NOUT, LERR, OK );

         // CSYTF2_RK

         SRNAMT = 'CSYTF2_RK';
         INFOT = 1;
         csytf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('CSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('CSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('CSYTF2_RK', INFOT, NOUT, LERR, OK );

         // CSYTRI_3

         SRNAMT = 'CSYTRI_3';
         INFOT = 1;
         csytri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('CSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('CSYTRI_3', INFOT, NOUT, LERR, OK );

         // CSYTRI_3X

         SRNAMT = 'CSYTRI_3X';
         INFOT = 1;
         csytri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('CSYTRI_3X', INFOT, NOUT, LERR, OK );

         // CSYTRS_3

         SRNAMT = 'CSYTRS_3';
         INFOT = 1;
         csytrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('CSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('CSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csytrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('CSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csytrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('CSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         csytrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('CSYTRS_3', INFOT, NOUT, LERR, OK );

         // CSYCON_3

         SRNAMT = 'CSYCON_3';
         INFOT = 1;
         csycon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csycon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csycon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('CSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, INFO);
         chkxer('CSYCON_3', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite packed matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // CSPTRF

         SRNAMT = 'CSPTRF';
         INFOT = 1;
         csptrf('/', 0, A, IP, INFO );
         chkxer('CSPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csptrf('U', -1, A, IP, INFO );
         chkxer('CSPTRF', INFOT, NOUT, LERR, OK );

         // CSPTRI

         SRNAMT = 'CSPTRI';
         INFOT = 1;
         csptri('/', 0, A, IP, W, INFO );
         chkxer('CSPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csptri('U', -1, A, IP, W, INFO );
         chkxer('CSPTRI', INFOT, NOUT, LERR, OK );

         // CSPTRS

         SRNAMT = 'CSPTRS';
         INFOT = 1;
         csptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('CSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('CSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('CSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('CSPTRS', INFOT, NOUT, LERR, OK );

         // CSPRFS

         SRNAMT = 'CSPRFS';
         INFOT = 1;
         csprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('CSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('CSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('CSPRFS', INFOT, NOUT, LERR, OK );

         // CSPCON

         SRNAMT = 'CSPCON';
         INFOT = 1;
         cspcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('CSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         cspcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('CSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         cspcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('CSPCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm

         // CSYTRF_AA

         SRNAMT = 'CSYTRF_AA';
         INFOT = 1;
         csytrf_aa('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrf_aa('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('CSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytrf_aa('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('CSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf_aa('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('CSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrf_aa('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('CSYTRF_AA', INFOT, NOUT, LERR, OK );

         // CSYTRS_AA

         SRNAMT = 'CSYTRS_AA';
         INFOT = 1;
         csytrs_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrs_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csytrs_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csytrs_aa('U', 2, 1, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         csytrs_aa('U', 2, 1, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csytrs_aa('U', 0, 1, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csytrs_aa('U', 0, 1, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('CSYTRS_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // CSYTRF_AA_2STAGE

         SRNAMT = 'CSYTRF_AA_2STAGE';
         INFOT = 1;
         csytrf_aa_2stage('/', 0, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('CSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrf_aa_2stage('U', -1, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('CSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         csytrf_aa_2stage('U', 2, A, 1, A, 2, IP, IP, W, 1, INFO );
         chkxer('CSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         csytrf_aa_2stage('U', 2, A, 2, A, 1, IP, IP, W, 1, INFO );
         chkxer('CSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         csytrf_aa_2stage('U', 2, A, 2, A, 8, IP, IP, W, 0, INFO );
         chkxer('CSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );

         // CHETRS_AA_2STAGE

         SRNAMT = 'CSYTRS_AA_2STAGE';
         INFOT = 1;
         csytrs_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         csytrs_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         csytrs_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         csytrs_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         csytrs_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         csytrs_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, INFO );
         chkxer('CSYTRS_AA_STAGE', INFOT, NOUT, LERR, OK );

      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of CERRSY

      }
