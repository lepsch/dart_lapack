      SUBROUTINE ZERRHE( PATH, NUNIT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             PATH;
      int                NUNIT;
      // ..

*  =====================================================================


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
      COMPLEX*16         A( NMAX, NMAX ), AF( NMAX, NMAX ), B( NMAX ), E( NMAX ), W( 2*NMAX ), X( NMAX )
      // ..
      // .. External Functions ..
      bool               LSAMEN;
      // EXTERNAL LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAESM, CHKXER, ZHECON, ZHECON_3, ZHECON_ROOK, ZHERFS, ZHETF2, ZHETF2_RK, ZHETF2_ROOK, ZHETRF, ZHETRF_RK, ZHETRF_ROOK, ZHETRF_AA, ZHETRF_AA_2STAGE, ZHETRI, ZHETRI_3, ZHETRI_3X, ZHETRI_ROOK, ZHETRI2, ZHETRI2X, ZHETRS, ZHETRS_3, ZHETRS_ROOK, ZHETRS_AA, ZHETRS_AA_2STAGE, ZHPCON, ZHPRFS, ZHPTRF, ZHPTRI, ZHPTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NOUT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX
      // ..
      // .. Executable Statements ..

      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
      C2 = PATH( 2: 3 )

      // Set the variables to innocuous values.

      for (J = 1; J <= NMAX; J++) { // 20
         for (I = 1; I <= NMAX; I++) { // 10
            A( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )             AF( I, J ) = DCMPLX( 1.D0 / DBLE( I+J ), -1.D0 / DBLE( I+J ) )
         } // 10
         B( J ) = 0.D0
         E( J ) = 0.D0
         R1( J ) = 0.D0
         R2( J ) = 0.D0
         W( J ) = 0.D0
         X( J ) = 0.D0
         IP( J ) = J
      } // 20
      ANRM = 1.0D0
      OK = .TRUE.

      if ( LSAMEN( 2, C2, 'HE' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // ZHETRF

         SRNAMT = 'ZHETRF'
         INFOT = 1
         zhetrf('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrf('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetrf('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF', INFOT, NOUT, LERR, OK );

         // ZHETF2

         SRNAMT = 'ZHETF2'
         INFOT = 1
         zhetf2('/', 0, A, 1, IP, INFO );
         chkxer('ZHETF2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetf2('U', -1, A, 1, IP, INFO );
         chkxer('ZHETF2', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetf2('U', 2, A, 1, IP, INFO );
         chkxer('ZHETF2', INFOT, NOUT, LERR, OK );

         // ZHETRI

         SRNAMT = 'ZHETRI'
         INFOT = 1
         zhetri('/', 0, A, 1, IP, W, INFO );
         chkxer('ZHETRI', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri('U', -1, A, 1, IP, W, INFO );
         chkxer('ZHETRI', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri('U', 2, A, 1, IP, W, INFO );
         chkxer('ZHETRI', INFOT, NOUT, LERR, OK );

         // ZHETRI2

         SRNAMT = 'ZHETRI2'
         INFOT = 1
         zhetri2('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri2('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri2('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2', INFOT, NOUT, LERR, OK );

         // ZHETRI2X

         SRNAMT = 'ZHETRI2X'
         INFOT = 1
         zhetri2x('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri2x('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri2x('U', 2, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRI2X', INFOT, NOUT, LERR, OK );

         // ZHETRS

         SRNAMT = 'ZHETRS'
         INFOT = 1
         zhetrs('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrs('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhetrs('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhetrs('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZHETRS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetrs('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZHETRS', INFOT, NOUT, LERR, OK );

         // ZHERFS

         SRNAMT = 'ZHERFS'
         INFOT = 1
         zherfs('/', 0, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zherfs('U', -1, 0, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zherfs('U', 0, -1, A, 1, AF, 1, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zherfs('U', 2, 1, A, 1, AF, 2, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zherfs('U', 2, 1, A, 2, AF, 1, IP, B, 2, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zherfs('U', 2, 1, A, 2, AF, 2, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );
         INFOT = 12
         zherfs('U', 2, 1, A, 2, AF, 2, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHERFS', INFOT, NOUT, LERR, OK );

         // ZHECON

         SRNAMT = 'ZHECON'
         INFOT = 1
         zhecon('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhecon('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhecon('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhecon('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHECON', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HR' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with rook
         // (bounded Bunch-Kaufman) diagonal pivoting method.

         // ZHETRF_ROOK

         SRNAMT = 'ZHETRF_ROOK'
         INFOT = 1
         zhetrf_rook('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrf_rook('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetrf_rook('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf_rook('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf_rook('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF_ROOK', INFOT, NOUT, LERR, OK );

         // ZHETF2_ROOK

         SRNAMT = 'ZHETF2_ROOK'
         INFOT = 1
         zhetf2_rook('/', 0, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetf2_rook('U', -1, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetf2_rook('U', 2, A, 1, IP, INFO );
         chkxer('ZHETF2_ROOK', INFOT, NOUT, LERR, OK );

         // ZHETRI_ROOK

         SRNAMT = 'ZHETRI_ROOK'
         INFOT = 1
         zhetri_rook('/', 0, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri_rook('U', -1, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri_rook('U', 2, A, 1, IP, W, INFO );
         chkxer('ZHETRI_ROOK', INFOT, NOUT, LERR, OK );

         // ZHETRS_ROOK

         SRNAMT = 'ZHETRS_ROOK'
         INFOT = 1
         zhetrs_rook('/', 0, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrs_rook('U', -1, 0, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhetrs_rook('U', 0, -1, A, 1, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhetrs_rook('U', 2, 1, A, 1, IP, B, 2, INFO );
         chkxer('ZHETRS_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetrs_rook('U', 2, 1, A, 2, IP, B, 1, INFO );
         chkxer('ZHETRS_ROOK', INFOT, NOUT, LERR, OK );

         // ZHECON_ROOK

         SRNAMT = 'ZHECON_ROOK'
         INFOT = 1
         zhecon_rook('/', 0, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhecon_rook('U', -1, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhecon_rook('U', 2, A, 1, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhecon_rook('U', 1, A, 1, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHECON_ROOK', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HK' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with rook
         // (bounded Bunch-Kaufman) pivoting with the new storage
         // format for factors L ( or U) and D.

         // L (or U) is stored in A, diagonal of D is stored on the
         // diagonal of A, subdiagonal of D is stored in a separate array E.

         // ZHETRF_RK

         SRNAMT = 'ZHETRF_RK'
         INFOT = 1
         zhetrf_rk('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrf_rk('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetrf_rk('U', 2, A, 1, E, IP, W, 4, INFO );
         chkxer('ZHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetrf_rk('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZHETRF_RK', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetrf_rk('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZHETRF_RK', INFOT, NOUT, LERR, OK );

         // ZHETF2_RK

         SRNAMT = 'ZHETF2_RK'
         INFOT = 1
         zhetf2_rk('/', 0, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetf2_rk('U', -1, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetf2_rk('U', 2, A, 1, E, IP, INFO );
         chkxer('ZHETF2_RK', INFOT, NOUT, LERR, OK );

         // ZHETRI_3

         SRNAMT = 'ZHETRI_3'
         INFOT = 1
         zhetri_3('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri_3('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri_3('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetri_3('U', 0, A, 1, E, IP, W, 0, INFO );
         chkxer('ZHETRI_3', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetri_3('U', 0, A, 1, E, IP, W, -2, INFO );
         chkxer('ZHETRI_3', INFOT, NOUT, LERR, OK );

         // ZHETRI_3X

         SRNAMT = 'ZHETRI_3X'
         INFOT = 1
         zhetri_3x('/', 0, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetri_3x('U', -1, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetri_3x('U', 2, A, 1, E, IP, W, 1, INFO );
         chkxer('ZHETRI_3X', INFOT, NOUT, LERR, OK );

         // ZHETRS_3

         SRNAMT = 'ZHETRS_3'
         INFOT = 1
         zhetrs_3('/', 0, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrs_3('U', -1, 0, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhetrs_3('U', 0, -1, A, 1, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhetrs_3('U', 2, 1, A, 1, E, IP, B, 2, INFO );
         chkxer('ZHETRS_3', INFOT, NOUT, LERR, OK );
         INFOT = 9
         zhetrs_3('U', 2, 1, A, 2, E, IP, B, 1, INFO );
         chkxer('ZHETRS_3', INFOT, NOUT, LERR, OK );

         // ZHECON_3

         SRNAMT = 'ZHECON_3'
         INFOT = 1
         zhecon_3('/', 0, A, 1,  E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhecon_3('U', -1, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhecon_3('U', 2, A, 1, E, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHECON_3', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhecon_3('U', 1, A, 1, E, IP, -1.0D0, RCOND, W, INFO);
         chkxer('ZHECON_3', INFOT, NOUT, LERR, OK );

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite matrix with Aasen's algorithm.

      } else if ( LSAMEN( 2, C2, 'HA' ) ) {

         // ZHETRF_AA

         SRNAMT = 'ZHETRF_AA'
         INFOT = 1
         zhetrf_aa('/', 0, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrf_aa('U', -1, A, 1, IP, W, 1, INFO );
         chkxer('ZHETRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetrf_aa('U', 2, A, 1, IP, W, 4, INFO );
         chkxer('ZHETRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf_aa('U', 0, A, 1, IP, W, 0, INFO );
         chkxer('ZHETRF_AA', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrf_aa('U', 0, A, 1, IP, W, -2, INFO );
         chkxer('ZHETRF_AA', INFOT, NOUT, LERR, OK );

         // ZHETRS_AA

         SRNAMT = 'ZHETRS_AA'
         INFOT = 1
         zhetrs_aa('/', 0, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrs_aa('U', -1, 0, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhetrs_aa('U', 0, -1, A, 1, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhetrs_aa('U', 2, 1, A, 1, IP, B, 2, W, 1, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhetrs_aa('U', 2, 1, A, 2, IP, B, 1, W, 1, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhetrs_aa('U', 0, 1, A, 1, IP, B, 1, W, 0, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhetrs_aa('U', 0, 1, A, 1, IP, B, 1, W, -2, INFO );
         chkxer('ZHETRS_AA', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'S2' ) ) {

         // Test error exits of the routines that use factorization
         // of a symmetric indefinite matrix with Aasen's algorithm.

         // ZHETRF_AA_2STAGE

         SRNAMT = 'ZHETRF_AA_2STAGE'
         INFOT = 1
         zhetrf_aa_2stage('/', 0, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrf_aa_2stage('U', -1, A, 1, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 4
         zhetrf_aa_2stage('U', 2, A, 1, A, 2, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 6
         zhetrf_aa_2stage('U', 2, A, 2, A, 1, IP, IP, W, 1, INFO );
         chkxer('ZHETRF_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhetrf_aa_2stage('U', 2, A, 2, A, 8, IP, IP, W, 0, INFO );
         chkxer('ZHETRF_AA_2STAGE', INFOT, NOUT, LERR, OK );

         // ZHETRS_AA_2STAGE

         SRNAMT = 'ZHETRS_AA_2STAGE'
         INFOT = 1
         zhetrs_aa_2stage('/', 0, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhetrs_aa_2stage('U', -1, 0, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhetrs_aa_2stage('U', 0, -1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhetrs_aa_2stage('U', 2, 1, A, 1, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhetrs_aa_2stage('U', 2, 1, A, 2, A, 1, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_2STAGE', INFOT, NOUT, LERR, OK );
         INFOT = 11
         zhetrs_aa_2stage('U', 2, 1, A, 2, A, 8, IP, IP, B, 1, INFO );
         chkxer('ZHETRS_AA_STAGE', INFOT, NOUT, LERR, OK );

      } else if ( LSAMEN( 2, C2, 'HP' ) ) {

         // Test error exits of the routines that use factorization
         // of a Hermitian indefinite packed matrix with partial
         // (Bunch-Kaufman) diagonal pivoting method.

         // ZHPTRF

         SRNAMT = 'ZHPTRF'
         INFOT = 1
         zhptrf('/', 0, A, IP, INFO );
         chkxer('ZHPTRF', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhptrf('U', -1, A, IP, INFO );
         chkxer('ZHPTRF', INFOT, NOUT, LERR, OK );

         // ZHPTRI

         SRNAMT = 'ZHPTRI'
         INFOT = 1
         zhptri('/', 0, A, IP, W, INFO );
         chkxer('ZHPTRI', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhptri('U', -1, A, IP, W, INFO );
         chkxer('ZHPTRI', INFOT, NOUT, LERR, OK );

         // ZHPTRS

         SRNAMT = 'ZHPTRS'
         INFOT = 1
         zhptrs('/', 0, 0, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhptrs('U', -1, 0, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhptrs('U', 0, -1, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', INFOT, NOUT, LERR, OK );
         INFOT = 7
         zhptrs('U', 2, 1, A, IP, B, 1, INFO );
         chkxer('ZHPTRS', INFOT, NOUT, LERR, OK );

         // ZHPRFS

         SRNAMT = 'ZHPRFS'
         INFOT = 1
         zhprfs('/', 0, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhprfs('U', -1, 0, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 3
         zhprfs('U', 0, -1, A, AF, IP, B, 1, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 8
         zhprfs('U', 2, 1, A, AF, IP, B, 1, X, 2, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', INFOT, NOUT, LERR, OK );
         INFOT = 10
         zhprfs('U', 2, 1, A, AF, IP, B, 2, X, 1, R1, R2, W, R, INFO );
         chkxer('ZHPRFS', INFOT, NOUT, LERR, OK );

         // ZHPCON

         SRNAMT = 'ZHPCON'
         INFOT = 1
         zhpcon('/', 0, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', INFOT, NOUT, LERR, OK );
         INFOT = 2
         zhpcon('U', -1, A, IP, ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', INFOT, NOUT, LERR, OK );
         INFOT = 5
         zhpcon('U', 1, A, IP, -ANRM, RCOND, W, INFO );
         chkxer('ZHPCON', INFOT, NOUT, LERR, OK );
      }

      // Print a summary line.

      alaesm(PATH, OK, NOUT );

      RETURN

      // End of ZERRHE

      }
