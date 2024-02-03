      SUBROUTINE ZLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B, LDB, ISEED, INFO )

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
      COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      const              ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) ;
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
      // EXTERNAL XERBLA, ZGBMV, ZGEMM, ZHBMV, ZHEMM, ZHPMV, ZLACPY, ZLARNV, ZSBMV, ZSPMV, ZSYMM, ZTBMV, ZTPMV, ZTRMM
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
      if ( .NOT.LSAME( C1, 'Zomplex precision' ) ) {
         INFO = -1
      } else if ( .NOT.( LSAME( XTYPE, 'N' ) .OR. LSAME( XTYPE, 'C' ) ) ) {
         INFO = -2
      } else if ( ( SYM .OR. TRI ) .AND. .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( ( GEN .OR. QRS ) .AND. .NOT. ( TRAN .OR. LSAME( TRANS, 'N' ) ) ) {
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
         CALL XERBLA( 'ZLARHS', -INFO )
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
         DO 10 J = 1, NRHS
            CALL ZLARNV( 2, ISEED, N, X( 1, J ) )
   10    CONTINUE
      }

      // Multiply X by op(A) using an appropriate
      // matrix multiply routine.

      if ( LSAMEN( 2, C2, 'GE' ) .OR. LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C2, 'RQ' ) ) {

         // General matrix

         CALL ZGEMM( TRANS, 'N', MB, NRHS, NX, ONE, A, LDA, X, LDX, ZERO, B, LDB )

      } else if ( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'HE' ) ) {

         // Hermitian matrix, 2-D storage

         CALL ZHEMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB )

      } else if ( LSAMEN( 2, C2, 'SY' ) ) {

         // Symmetric matrix, 2-D storage

         CALL ZSYMM( 'Left', UPLO, N, NRHS, ONE, A, LDA, X, LDX, ZERO, B, LDB )

      } else if ( LSAMEN( 2, C2, 'GB' ) ) {

         // General matrix, band storage

         DO 20 J = 1, NRHS
            CALL ZGBMV( TRANS, M, N, KL, KU, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   20    CONTINUE

      } else if ( LSAMEN( 2, C2, 'PB' ) .OR. LSAMEN( 2, C2, 'HB' ) ) {

         // Hermitian matrix, band storage

         DO 30 J = 1, NRHS
            CALL ZHBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   30    CONTINUE

      } else if ( LSAMEN( 2, C2, 'SB' ) ) {

         // Symmetric matrix, band storage

         DO 40 J = 1, NRHS
            CALL ZSBMV( UPLO, N, KL, ONE, A, LDA, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   40    CONTINUE

      } else if ( LSAMEN( 2, C2, 'PP' ) .OR. LSAMEN( 2, C2, 'HP' ) ) {

         // Hermitian matrix, packed storage

         DO 50 J = 1, NRHS
            CALL ZHPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   50    CONTINUE

      } else if ( LSAMEN( 2, C2, 'SP' ) ) {

         // Symmetric matrix, packed storage

         DO 60 J = 1, NRHS
            CALL ZSPMV( UPLO, N, ONE, A, X( 1, J ), 1, ZERO, B( 1, J ), 1 )
   60    CONTINUE

      } else if ( LSAMEN( 2, C2, 'TR' ) ) {

         // Triangular matrix.  Note that for triangular matrices,
            // KU = 1 => non-unit triangular
            // KU = 2 => unit triangular

         CALL ZLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         if ( KU.EQ.2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         CALL ZTRMM( 'Left', UPLO, TRANS, DIAG, N, NRHS, ONE, A, LDA, B, LDB )

      } else if ( LSAMEN( 2, C2, 'TP' ) ) {

         // Triangular matrix, packed storage

         CALL ZLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         if ( KU.EQ.2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         DO 70 J = 1, NRHS
            CALL ZTPMV( UPLO, TRANS, DIAG, N, A, B( 1, J ), 1 )
   70    CONTINUE

      } else if ( LSAMEN( 2, C2, 'TB' ) ) {

         // Triangular matrix, banded storage

         CALL ZLACPY( 'Full', N, NRHS, X, LDX, B, LDB )
         if ( KU.EQ.2 ) {
            DIAG = 'U'
         } else {
            DIAG = 'N'
         }
         DO 80 J = 1, NRHS
            CALL ZTBMV( UPLO, TRANS, DIAG, N, KL, A, LDA, B( 1, J ), 1 )
   80    CONTINUE

      } else {

         // If none of the above, set INFO = -1 and return

         INFO = -1
         CALL XERBLA( 'ZLARHS', -INFO )
      }

      RETURN

      // End of ZLARHS

      }
