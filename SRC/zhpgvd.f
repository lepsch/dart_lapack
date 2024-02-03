      SUBROUTINE ZHPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * ), W( * );
      COMPLEX*16         AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                J, LIWMIN, LRWMIN, LWMIN, NEIG;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHPEVD, ZHPGST, ZPPTRF, ZTPMV, ZTPSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      INFO = 0
      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
         INFO = -9
      }

      if ( INFO.EQ.0 ) {
         if ( N.LE.1 ) {
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         } else {
            if ( WANTZ ) {
               LWMIN = 2*N
               LRWMIN = 1 + 5*N + 2*N**2
               LIWMIN = 3 + 5*N
            } else {
               LWMIN = N
               LRWMIN = N
               LIWMIN = 1
            }
         }

         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         } else if ( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) {
            INFO = -13
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -15
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZHPGVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a Cholesky factorization of B.

      zpptrf(UPLO, N, BP, INFO );
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      zhpgst(ITYPE, UPLO, N, AP, BP, INFO );
      zhpevd(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO );
      LWMIN = INT( MAX( DBLE( LWMIN ), DBLE( WORK( 1 ) ) ) )
      LRWMIN = INT( MAX( DBLE( LRWMIN ), DBLE( RWORK( 1 ) ) ) )
      LIWMIN = INT( MAX( DBLE( LIWMIN ), DBLE( IWORK( 1 ) ) ) )

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         IF( INFO.GT.0 ) NEIG = INFO - 1
         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            DO 10 J = 1, NEIG
               ztpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
   10       CONTINUE

         } else if ( ITYPE.EQ.3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H *y

            if ( UPPER ) {
               TRANS = 'C'
            } else {
               TRANS = 'N'
            }

            DO 20 J = 1, NEIG
               ztpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
   20       CONTINUE
         }
      }

      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
      RETURN

      // End of ZHPGVD

      }
