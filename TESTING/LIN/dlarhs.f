      SUBROUTINE DLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, ISEED, INFO );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS, UPLO, XTYPE;
      String             PATH;
      int                INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      double             A( LDA, * ), B( LDB, * ), X( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               BAND, GEN, NOTRAN, QRS, SYM, TRAN, TRI;
      String             C1, DIAG;
      String             C2;
      int                J, MB, NX;
      // ..
      // .. External Functions ..
      bool               LSAME, LSAMEN;
      // EXTERNAL LSAME, LSAMEN
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGBMV, DGEMM, DLACPY, DLARNV, DSBMV, DSPMV, DSYMM, DTBMV, DTPMV, DTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      C1 = PATH( 1: 1 );
      C2 = PATH( 2: 3 );
      TRAN = LSAME( TRANS, 'T' ) || LSAME( TRANS, 'C' );
      NOTRAN = !TRAN;
      GEN = LSAME( PATH( 2: 2 ), 'G' );
      QRS = LSAME( PATH( 2: 2 ), 'Q' ) || LSAME( PATH( 3: 3 ), 'Q' );
      SYM = LSAME( PATH( 2: 2 ), 'P' ) || LSAME( PATH( 2: 2 ), 'S' );
      TRI = LSAME( PATH( 2: 2 ), 'T' );
      BAND = LSAME( PATH( 3: 3 ), 'B' );
      if ( !LSAME( C1, 'double;
         ' ) ) {
         INFO = -1;
      } else if ( !( LSAME( XTYPE, 'N' ) || LSAME( XTYPE, 'C' ) ) ) {
         INFO = -2;
      } else if ( ( SYM || TRI ) && !( LSAME( UPLO, 'U' ) || LSAME( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( ( GEN || QRS ) && !( TRAN || LSAME( TRANS, 'N' ) ) ) {
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
      } else if ( ( !BAND && LDA < MAX( 1, M ) ) || ( BAND && ( SYM || TRI ) && LDA < KL+1 ) || ( BAND && GEN && LDA < KL+KU+1 ) ) {
         INFO = -11;
      } else if ( ( NOTRAN && LDX < MAX( 1, N ) ) || ( TRAN && LDX < MAX( 1, M ) ) ) {
         INFO = -13;
      } else if ( ( NOTRAN && LDB < MAX( 1, M ) ) || ( TRAN && LDB < MAX( 1, N ) ) ) {
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
      if ( !LSAME( XTYPE, 'C' ) ) {
         for (J = 1; J <= NRHS; J++) { // 10
            dlarnv(2, ISEED, N, X( 1, J ) );
         } // 10
      }

      // Multiply X by op(A) using an appropriate
      // matrix multiply routine.

      if ( LSAMEN( 2, C2, 'GE' ) || LSAMEN( 2, C2, 'QR' ) || LSAMEN( 2, C2, 'LQ' ) || LSAMEN( 2, C2, 'QL' ) || LSAMEN( 2, C2, 'RQ' ) ) {

         // General matrix

         dgemm(TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( LSAMEN( 2, C2, 'PO' ) || LSAMEN( 2, C2, 'SY' ) ) {

         // Symmetric matrix, 2-D storage

         dsymm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // General matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 20
            dgbmv(TRANS, MB, NX, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 20

      } else if ( LSAMEN( 2, C2, 'PB' ) ) {

         // Symmetric matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 30
            dsbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 30

      } else if ( LSAMEN( 2, C2, 'PP' ) || LSAMEN( 2, C2, 'SP' ) ) {

         // Symmetric matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 40
            dspmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 40

      } else if ( LSAMEN( 2, C2, 'TR' ) ) {

         // Triangular matrix.  Note that for triangular matrices,
            // KU = 1 => non-unit triangular
            // KU = 2 => unit triangular

         dlacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U';
         } else {
            DIAG = 'N';
         }
         dtrmm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      } else if ( LSAMEN( 2, C2, 'TP' ) ) {

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

      } else if ( LSAMEN( 2, C2, 'TB' ) ) {

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

      return;

      // End of DLARHS

      }
