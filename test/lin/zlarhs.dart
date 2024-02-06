      void zlarhs(PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, ISEED, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS, UPLO, XTYPE;
      String             PATH;
      int                INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS;
      int                ISEED( 4 );
      Complex         A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      Complex         ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      bool               BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI;
      String             C1, DIAG;
      String             C2;
      int                J, MB, NX;
      // ..
      // .. External Functions ..
      //- bool               lsame, LSAMEN;
      // EXTERNAL lsame, LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGBMV, ZGEMM, ZHBMV, ZHEMM, ZHPMV, ZLACPY, ZLARNV, ZSBMV, ZSPMV, ZSYMM, ZTBMV, ZTPMV, ZTRMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX

      // Test the input parameters.

      INFO = 0;
      C1 = PATH( 1: 1 );
      C2 = PATH( 2: 3 );
      TRAN = lsame( TRANS, 'T' ) || lsame( TRANS, 'C' );
      NOTRAN = !TRAN;
      GEN = lsame( PATH( 2: 2 ), 'G' );
      QRS = lsame( PATH( 2: 2 ), 'Q' ) || lsame( PATH( 3: 3 ), 'Q' );
      SYM = lsame( PATH( 2: 2 ), 'P' ) || lsame( PATH( 2: 2 ), 'S' ) || lsame( PATH( 2: 2 ), 'H' );
      TRI = lsame( PATH( 2: 2 ), 'T' );
      BAND = lsame( PATH( 3: 3 ), 'B' );
      if ( !lsame( C1, 'Zomplex precision' ) ) {
         INFO = -1;
      } else if ( !( lsame( XTYPE, 'N' ) || lsame( XTYPE, 'C' ) ) ) {
         INFO = -2;
      } else if ( ( SYM || TRI ) && !( lsame( UPLO, 'U' ) || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( ( GEN || QRS ) && !( TRAN || lsame( TRANS, 'N' ) ) ) {
         INFO = -4;
      } else if ( M < 0 ) {
         INFO = -5;
      } else if ( N < 0 ) {
         INFO = -6;
      } else if ( BAND && KL < 0 ) {
         INFO = -7;
      } else if ( BAND && KU < 0 ) {
         INFO = -8;
      } else if ( NRHS < 0 ) {
         INFO = -9;
      } else if ( ( !BAND && LDA < max( 1, M ) ) || ( BAND && ( SYM || TRI ) && LDA < KL+1 ) || ( BAND && GEN && LDA < KL+KU+1 ) ) {
         INFO = -11;
      } else if ( ( NOTRAN && LDX < max( 1, N ) ) || ( TRAN && LDX < max( 1, M ) ) ) {
         INFO = -13;
      } else if ( ( NOTRAN && LDB < max( 1, M ) ) || ( TRAN && LDB < max( 1, N ) ) ) {
         INFO = -15;
      }
      if ( INFO != 0 ) {
         xerbla('ZLARHS', -INFO );
         return;
      }

      // Initialize X to NRHS random vectors unless XTYPE = 'C'.

      if ( TRAN ) {
         NX = M;
         MB = N;
      } else {
         NX = N;
         MB = M;
      }
      if ( !lsame( XTYPE, 'C' ) ) {
         for (J = 1; J <= NRHS; J++) { // 10
            zlarnv(2, ISEED, N, X( 1, J ) );
         } // 10
      }

      // Multiply X by op(A) using an appropriate
      // matrix multiply routine.

      if ( lsamen( 2, C2, 'GE' ) || lsamen( 2, C2, 'QR' ) || lsamen( 2, C2, 'LQ' ) || lsamen( 2, C2, 'QL' ) || lsamen( 2, C2, 'RQ' ) ) {

         // General matrix

         zgemm(TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( lsamen( 2, C2, 'PO' ) || lsamen( 2, C2, 'HE' ) ) {

         // Hermitian matrix, 2-D storage

         zhemm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( lsamen( 2, C2, 'SY' ) ) {

         // Symmetric matrix, 2-D storage

         zsymm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // General matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 20
            zgbmv(TRANS, M, N, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 20

      } else if ( lsamen( 2, C2, 'PB' ) || lsamen( 2, C2, 'HB' ) ) {

         // Hermitian matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 30
            zhbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 30

      } else if ( lsamen( 2, C2, 'SB' ) ) {

         // Symmetric matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 40
            zsbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 40

      } else if ( lsamen( 2, C2, 'PP' ) || lsamen( 2, C2, 'HP' ) ) {

         // Hermitian matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 50
            zhpmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 50

      } else if ( lsamen( 2, C2, 'SP' ) ) {

         // Symmetric matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 60
            zspmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 60

      } else if ( lsamen( 2, C2, 'TR' ) ) {

         // Triangular matrix.  Note that for triangular matrices,
            // KU = 1 => non-unit triangular
            // KU = 2 => unit triangular

         zlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         ztrmm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // Triangular matrix, packed storage

         zlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         for (J = 1; J <= NRHS; J++) { // 70
            ztpmv(UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 );
         } // 70

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // Triangular matrix, banded storage

         zlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         for (J = 1; J <= NRHS; J++) { // 80
            ztbmv(UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 );
         } // 80

      } else {

         // If none of the above, set INFO = -1 and return;

         INFO = -1;
         xerbla('ZLARHS', -INFO );
      }

      return;
      }
