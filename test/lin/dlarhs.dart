      void dlarhs(final int PATH, final int XTYPE, final int UPLO, final int TRANS, final int M, final int N, final int KL, final int KU, final int NRHS, final Matrix<double> A, final int LDA, final Matrix<double> X, final int LDX, final Matrix<double> B, final int LDB, final Array<int> ISEED, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TRANS, UPLO, XTYPE;
      String             PATH;
      int                INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS;
      int                ISEED( 4 );
      double             A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
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
      // EXTERNAL DGBMV, DGEMM, DLACPY, DLARNV, DSBMV, DSPMV, DSYMM, DTBMV, DTPMV, DTRMM, XERBLA
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
      SYM = lsame( PATH( 2: 2 ), 'P' ) || lsame( PATH( 2: 2 ), 'S' );
      TRI = lsame( PATH( 2: 2 ), 'T' );
      BAND = lsame( PATH( 3: 3 ), 'B' );
      if ( !lsame( C1, 'double;
         ' ) ) {
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
         xerbla('DLARHS', -INFO );
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
            dlarnv(2, ISEED, N, X( 1, J ) );
         } // 10
      }

      // Multiply X by op(A) using an appropriate
      // matrix multiply routine.

      if ( lsamen( 2, C2, 'GE' ) || lsamen( 2, C2, 'QR' ) || lsamen( 2, C2, 'LQ' ) || lsamen( 2, C2, 'QL' ) || lsamen( 2, C2, 'RQ' ) ) {

         // General matrix

         dgemm(TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( lsamen( 2, C2, 'PO' ) || lsamen( 2, C2, 'SY' ) ) {

         // Symmetric matrix, 2-D storage

         dsymm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( lsamen( 2, C2, 'GB' ) ) {

         // General matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 20
            dgbmv(TRANS, MB, NX, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 20

      } else if ( lsamen( 2, C2, 'PB' ) ) {

         // Symmetric matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 30
            dsbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 30

      } else if ( lsamen( 2, C2, 'PP' ) || lsamen( 2, C2, 'SP' ) ) {

         // Symmetric matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 40
            dspmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 40

      } else if ( lsamen( 2, C2, 'TR' ) ) {

         // Triangular matrix.  Note that for triangular matrices,
         //    KU = 1 => non-unit triangular
         //    KU = 2 => unit triangular

         dlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         dtrmm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      } else if ( lsamen( 2, C2, 'TP' ) ) {

         // Triangular matrix, packed storage

         dlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         for (J = 1; J <= NRHS; J++) { // 50
            dtpmv(UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 );
         } // 50

      } else if ( lsamen( 2, C2, 'TB' ) ) {

         // Triangular matrix, banded storage

         dlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         for (J = 1; J <= NRHS; J++) { // 60
            dtbmv(UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 );
         } // 60

      } else {

         // If PATH is none of the above, return with an error code.

         INFO = -1;
         xerbla('DLARHS', -INFO );
      }

      }
