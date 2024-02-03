      SUBROUTINE SSPGVD( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               AP( * ), BP( * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                J, LIWMIN, LWMIN, NEIG;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPPTRF, SSPEVD, SSPGST, STPMV, STPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

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
            LIWMIN = 1
            LWMIN = 1
         } else {
            if ( WANTZ ) {
               LIWMIN = 3 + 5*N
               LWMIN = 1 + 6*N + 2*N**2
            } else {
               LIWMIN = 1
               LWMIN = 2*N
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN
         if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SSPGVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form a Cholesky factorization of BP.

      spptrf(UPLO, N, BP, INFO );
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      sspgst(ITYPE, UPLO, N, AP, BP, INFO );
      sspevd(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO );
      LWMIN = INT( MAX( REAL( LWMIN ), REAL( WORK( 1 ) ) ) )
      LIWMIN = INT( MAX( REAL( LIWMIN ), REAL( IWORK( 1 ) ) ) )

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         IF( INFO.GT.0 ) NEIG = INFO - 1
         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T *y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            DO 10 J = 1, NEIG
               stpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
   10       CONTINUE

         } else if ( ITYPE.EQ.3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T *y

            if ( UPPER ) {
               TRANS = 'T'
            } else {
               TRANS = 'N'
            }

            DO 20 J = 1, NEIG
               stpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
   20       CONTINUE
         }
      }

      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SSPGVD

      }
