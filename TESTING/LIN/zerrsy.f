      void zerrsy(PATH, NUNIT ) {

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
      double             ANRM, RCOND;
      // ..
      // .. Local Arrays ..
      int                IP( NMAX );
      double             R( NMAX ), R1( NMAX ), R2( NMAX );
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX );
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZSPCON, ZSPRFS, ZSPTRF, ZSPTRI, ZSPTRS, ZSYCON, ZSYCON_3, ZSYCON_ROOK, ZSYRFS, ZSYTF2, ZSYTF2_RK, ZSYTF2_ROOK, ZSYTRF, ZSYTRF_RK, ZSYTRF_ROOK, ZSYTRI, ZSYTRI_3, ZSYTRI_3X, ZSYTRI_ROOK, ZSYTRI2, ZSYTRI2X, ZSYTRS, ZSYTRS_3, ZSYTRS_ROOK
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
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT;
      WRITE( NOUT, FMT = * );
      C2 = PATH( 2: 3 );

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.0 / DBLE( I+J ), -1.0 / DBLE( I+J ) );
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

         // ZSYTRF

         SRNAMT = 'ZSYTRF';
         INFOT = 1;
         zsytrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF', INFOT, NOUT, LERR, OK );

         // ZSYTF2

         SRNAMT = 'ZSYTF2';
         INFOT = 1;
         zsytf2('/', 0, A, 1, IP, INFO );
         chkxer('ZSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytf2('U', -1, A, 1, IP, INFO );
         chkxer('ZSYTF2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytf2('U', 2, A, 1, IP, INFO );
         chkxer('ZSYTF2', INFOT, NOUT, LERR, OK );

         // ZSYTRI

         SRNAMT = 'ZSYTRI';
         INFOT = 1;
         zsytri('/', 0, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri('U', -1, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri('U', 2, A, 1, IP, W, INFO );
         chkxer('ZSYTRI', INFOT, NOUT, LERR, OK );

         // ZSYTRI2

         SRNAMT = 'ZSYTRI2';
         INFOT = 1;
         zsytri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2', INFOT, NOUT, LERR, OK );

         // ZSYTRI2X

         SRNAMT = 'ZSYTRI2X';
         INFOT = 1;
         zsytri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRI2X', INFOT, NOUT, LERR, OK );

         // ZSYTRS

         SRNAMT = 'ZSYTRS';
         INFOT = 1;
         zsytrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsytrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsytrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZSYTRS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZSYTRS', INFOT, NOUT, LERR, OK );

         // ZSYRFS

         SRNAMT = 'ZSYRFS';
         INFOT = 1;
         zsyrfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsyrfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsyrfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsyrfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsyrfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );
         INFOT = 12;
         zsyrfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSYRFS', INFOT, NOUT, LERR, OK );

         // ZSYCON

         SRNAMT = 'ZSYCON';
         INFOT = 1;
         zsycon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsycon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsycon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zsycon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSYCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SR' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // ZSYTRF_ROOK

         SRNAMT = 'ZSYTRF_ROOK';
         INFOT = 1;
         zsytrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF_ROOK', INFOT, NOUT, LERR, OK );

         // ZSYTF2_ROOK

         SRNAMT = 'ZSYTF2_ROOK';
         INFOT = 1;
         zsytf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('ZSYTF2_ROOK', INFOT, NOUT, LERR, OK );

         // ZSYTRI_ROOK

         SRNAMT = 'ZSYTRI_ROOK';
         INFOT = 1;
         zsytri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('ZSYTRI_ROOK', INFOT, NOUT, LERR, OK );

         // ZSYTRS_ROOK

         SRNAMT = 'ZSYTRS_ROOK';
         INFOT = 1;
         zsytrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsytrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsytrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZSYTRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZSYTRS_ROOK', INFOT, NOUT, LERR, OK );

         // ZSYCON_ROOK

         SRNAMT = 'ZSYCON_ROOK';
         INFOT = 1;
         zsycon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsycon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsycon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zsycon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // ZSYTRF_RK

         SRNAMT = 'ZSYTRF_RK';
         INFOT = 1;
         zsytrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('ZSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZSYTRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZSYTRF_RK', INFOT, NOUT, LERR, OK );

         // ZSYTF2_RK

         SRNAMT = 'ZSYTF2_RK';
         INFOT = 1;
         zsytf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('ZSYTF2_RK', INFOT, NOUT, LERR, OK );

         // ZSYTRI_3

         SRNAMT = 'ZSYTRI_3';
         INFOT = 1;
         zsytri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZSYTRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZSYTRI_3', INFOT, NOUT, LERR, OK );

         // ZSYTRI_3X

         SRNAMT = 'ZSYTRI_3X';
         INFOT = 1;
         zsytri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZSYTRI_3X', INFOT, NOUT, LERR, OK );

         // ZSYTRS_3

         SRNAMT = 'ZSYTRS_3';
         INFOT = 1;
         zsytrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsytrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsytrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('ZSYTRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9;
         zsytrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('ZSYTRS_3', INFOT, NOUT, LERR, OK );

         // ZSYCON_3

         SRNAMT = 'ZSYCON_3';
         INFOT = 1;
         zsycon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsycon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsycon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSYCON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsycon_3('U', 1, A, 1, E, IP, -1.0, RCOND, W, INFO);
         chkxer('ZSYCON_3', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite packed matrix with partial
         // (Bunch-Kaufman) pivoting.

         // ZSPTRF

         SRNAMT = 'ZSPTRF';
         INFOT = 1;
         zsptrf('/', 0, A, IP, INFO );
         chkxer('ZSPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsptrf('U', -1, A, IP, INFO );
         chkxer('ZSPTRF', INFOT, NOUT, LERR, OK );

         // ZSPTRI

         SRNAMT = 'ZSPTRI';
         INFOT = 1;
         zsptri('/', 0, A, IP, W, INFO );
         chkxer('ZSPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsptri('U', -1, A, IP, W, INFO );
         chkxer('ZSPTRI', INFOT, NOUT, LERR, OK );

         // ZSPTRS

         SRNAMT = 'ZSPTRS';
         INFOT = 1;
         zsptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('ZSPTRS', INFOT, NOUT, LERR, OK );

         // ZSPRFS

         SRNAMT = 'ZSPRFS';
         INFOT = 1;
         zsprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zsprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZSPRFS', INFOT, NOUT, LERR, OK );

         // ZSPCON

         SRNAMT = 'ZSPCON';
         INFOT = 1;
         zspcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zspcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zspcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZSPCON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'SA' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZSYTRF_AA

         SRNAMT = 'ZSYTRF_AA';
         INFOT = 1;
         zsytrf_aa('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrf_aa('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytrf_aa('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf_aa('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZSYTRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrf_aa('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZSYTRF_AA', INFOT, NOUT, LERR, OK );

         // ZSYTRS_AA

         SRNAMT = 'ZSYTRS_AA';
         INFOT = 1;
         zsytrs_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrs_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsytrs_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsytrs_aa('U', 2, 1, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('ZSYTRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8;
         zsytrs_aa('U', 2, 1, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZSYTRS_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZSYTRF_AA_2STAGE

         SRNAMT = 'ZSYTRF_AA_2STAGE';
         INFOT = 1;
         zsytrf_aa_2stage('/', 0, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrf_aa_2stage('U', -1, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4;
         zsytrf_aa_2stage('U', 2, A, 1, A, 2, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6;
         zsytrf_aa_2stage('U', 2, A, 2, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10;
         zsytrf_aa_2stage('U', 2, A, 2, A, 8, IP, IP, W, 0, INFO );
         chkxer('ZSYTRF_AA_2STAGE', INFOT, NOUT, LERR, OK );

         // CHETRS_AA_2STAGE

         SRNAMT = 'ZSYTRS_AA_2STAGE';
         INFOT = 1;
         zsytrs_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2;
         zsytrs_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3;
         zsytrs_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5;
         zsytrs_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7;
         zsytrs_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11;
         zsytrs_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, INFO );
         chkxer('ZSYTRS_AA_STAGE', INFOT, NOUT, LERR, OK );

      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      return;

      // End of ZERRSY

      }
