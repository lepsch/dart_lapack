      SUBROUTINE CLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, ISEED, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS, UPLO, XTYPE;
      String             PATH;
      int                INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX            A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO
      const              ONE = ( 1.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ) ;
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
      // EXTERNAL CGBMV, CGEMM, CHBMV, CHEMM, CHPMV, CLACPY, CLARNV, CSBMV, CSPMV, CSYMM, CTBMV, CTPMV, CTRMM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      C1 = PATH( 1: 1 )
      C2 = PATH( 2: 3 )
      TRAN = LSAME( TRANS, 'T' ) .OR. LSAME( TRANS, 'C' )
      NOTRAN = .NOT.TRAN
      GEN = LSAME( PATH( 2: 2 ), 'G' )
      QRS = LSAME( PATH( 2: 2 ), 'Q' ) .OR. LSAME( PATH( 3: 3 ), 'Q' )
      SYM = LSAME( PATH( 2: 2 ), 'P' ) .OR. LSAME( PATH( 2: 2 ), 'S' ) .OR. LSAME( PATH( 2: 2 ), 'H' )
      TRI = LSAME( PATH( 2: 2 ), 'T' )
      BAND = LSAME( PATH( 3: 3 ), 'B' )
      if ( .NOT.LSAME( C1, 'Complex precision' ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( XTYPE, 'N' ) .OR. LSAME( XTYPE, 'C' ) ) ) {
         INFO = -2
      } else if ( ( SYM .OR. TRI ) .AND. .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( ( GEN.OR.QRS ) .AND. .NOT.( TRAN .OR. LSAME( TRANS, 'N' ) ) ) {
         INFO = -4
      } else if ( M.LT.0 ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( BAND .AND. KL.LT.0 ) {
         INFO = -7
      } else if ( BAND .AND. KU.LT.0 ) {
         INFO = -8
      } else if ( NRHS.LT.0 ) {
         INFO = -9
      } else if ( ( .NOT.BAND .AND. LDA.LT.MAX( 1, M ) ) .OR. ( BAND .AND. ( SYM .OR. TRI ) .AND. LDA.LT.KL+1 ) .OR. ( BAND .AND. GEN .AND. LDA.LT.KL+KU+1 ) ) {
         INFO = -11
      } else if ( ( NOTRAN .AND. LDX.LT.MAX( 1, N ) ) .OR. ( TRAN .AND. LDX.LT.MAX( 1, M ) ) ) {
         INFO = -13
      } else if ( ( NOTRAN .AND. LDB.LT.MAX( 1, M ) ) .OR. ( TRAN .AND. LDB.LT.MAX( 1, N ) ) ) {
         INFO = -15
      }
      if ( INFO.NE.0 ) {
         xerbla('CLARHS', -INFO );
         RETURN
      }

      // Initialize X to NRHS random vectors unless XTYPE = 'C'.

      if ( TRAN ) {
         NX = M
         MB = N
      } else {
         NX = N
         MB = M
      }
      if ( .NOT.LSAME( XTYPE, 'C' ) ) {
         for (J = 1; J <= NRHS; J++) { // 10
            clarnv(2, ISEED, N, X( 1, J ) );
         } // 10
      }

      // Multiply X by op(A) using an appropriate
      // matrix multiply routine.

      if ( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) {

         // General matrix

         cgemm(TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'HE' ) ) {

         // Hermitian matrix, 2-D storage

         chemm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // Symmetric matrix, 2-D storage

         csymm('Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB );

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // General matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 20
            cgbmv(TRANS, M, N, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 20

      } else if ( LSAMEN( 2, C2, 'PB' ) .OR. LSAMEN( 2, C2, 'HB' ) ) {

         // Hermitian matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 30
            chbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 30

      } else if ( LSAMEN( 2, C2, 'SB' ) ) {

         // Symmetric matrix, band storage

         for (J = 1; J <= NRHS; J++) { // 40
            csbmv(UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 40

      } else if ( LSAMEN( 2, C2, 'PP' ) .OR. LSAMEN( 2, C2, 'HP' ) ) {

         // Hermitian matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 50
            chpmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 50

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Symmetric matrix, packed storage

         for (J = 1; J <= NRHS; J++) { // 60
            cspmv(UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 );
         } // 60

      } else if ( LSAMEN( 2, C2, 'TR' ) ) {

         // Triangular matrix.  Note that for triangular matrices,
            // KU = 1 => non-unit triangular
            // KU = 2 => unit triangular

         clacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         ctrmm('Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB );

      } else if ( LSAMEN( 2, C2, 'TP' ) ) {

         // Triangular matrix, packed storage

         clacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         for (J = 1; J <= NRHS; J++) { // 70
            ctpmv(UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 );
         } // 70

      } else if ( LSAMEN( 2, C2, 'TB' ) ) {

         // Triangular matrix, banded storage

         clacpy('Full', N, NRHS, X, LDX, B, LDB );
         if ( KU == 2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         for (J = 1; J <= NRHS; J++) { // 80
            ctbmv(UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 );
         } // 80

      } else {

         // If none of the above, set INFO = -1 and return

         INFO = -1
         xerbla('CLARHS', -INFO );
      }

      RETURN

      // End of CLARHS

      }
