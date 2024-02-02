      SUBROUTINE DGEBAL( JOB, N, A, LDA, ILO, IHI, SCALE, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), SCALE( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   SCLFAC
      PARAMETER          ( SCLFAC = 2.0D+0 )
      DOUBLE PRECISION   FACTOR
      PARAMETER          ( FACTOR = 0.95D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOCONV, CANSWAP
      INTEGER            I, ICA, IRA, J, K, L
      DOUBLE PRECISION   C, CA, F, G, R, RA, S, SFMAX1, SFMAX2, SFMIN1,
     $                   SFMIN2
*     ..
*     .. External Functions ..
      LOGICAL            DISNAN, LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DNRM2
      EXTERNAL           DISNAN, LSAME, IDAMAX, DLAMCH, DNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEBAL', -INFO )
         RETURN
      END IF
*
*     Quick returns.
*
      IF( N.EQ.0 ) THEN
         ILO = 1
         IHI = 0
         RETURN
      END IF
*
      IF( LSAME( JOB, 'N' ) ) THEN
         DO I = 1, N
            SCALE( I ) = ONE
         END DO
         ILO = 1
         IHI = N
         RETURN
      END IF
*
*     Permutation to isolate eigenvalues if possible.
*
      K = 1
      L = N
*
      IF( .NOT.LSAME( JOB, 'S' ) ) THEN
*
*        Row and column exchange.
*
         NOCONV = .TRUE.
         DO WHILE( NOCONV )
*
*           Search for rows isolating an eigenvalue and push them down.
*
            NOCONV = .FALSE.
            DO I = L, 1, -1
               CANSWAP = .TRUE.
               DO J = 1, L
                  IF( I.NE.J .AND. A( I, J ).NE.ZERO ) THEN
                     CANSWAP = .FALSE.
                     EXIT
                  END IF
               END DO
*
               IF( CANSWAP ) THEN
                  SCALE( L ) = I
                  IF( I.NE.L ) THEN
                     CALL DSWAP( L, A( 1, I ), 1, A( 1, L ), 1 )
                     CALL DSWAP( N-K+1, A( I, K ), LDA, A( L, K ), LDA )
                  END IF
                  NOCONV = .TRUE.
*
                  IF( L.EQ.1 ) THEN
                     ILO = 1
                     IHI = 1
                     RETURN
                  END IF
*
                  L = L - 1
               END IF
            END DO
*
         END DO

         NOCONV = .TRUE.
         DO WHILE( NOCONV )
*
*           Search for columns isolating an eigenvalue and push them left.
*
            NOCONV = .FALSE.
            DO J = K, L
               CANSWAP = .TRUE.
               DO I = K, L
                  IF( I.NE.J .AND. A( I, J ).NE.ZERO ) THEN
                     CANSWAP = .FALSE.
                     EXIT
                  END IF
               END DO
*
               IF( CANSWAP ) THEN
                  SCALE( K ) = J
                  IF( J.NE.K ) THEN
                     CALL DSWAP( L, A( 1, J ), 1, A( 1, K ), 1 )
                     CALL DSWAP( N-K+1, A( J, K ), LDA, A( K, K ), LDA )
                  END IF
                  NOCONV = .TRUE.
*
                  K = K + 1
               END IF
            END DO
*
         END DO
*
      END IF
*
*     Initialize SCALE for non-permuted submatrix.
*
      DO I = K, L
         SCALE( I ) = ONE
      END DO
*
*     If we only had to permute, we are done.
*
      IF( LSAME( JOB, 'P' ) ) THEN
         ILO = K
         IHI = L
         RETURN
      END IF
*
*     Balance the submatrix in rows K to L.
*
*     Iterative loop for norm reduction.
*
      SFMIN1 = DLAMCH( 'S' ) / DLAMCH( 'P' )
      SFMAX1 = ONE / SFMIN1
      SFMIN2 = SFMIN1*SCLFAC
      SFMAX2 = ONE / SFMIN2
*
      NOCONV = .TRUE.
      DO WHILE( NOCONV )
         NOCONV = .FALSE.
*
         DO I = K, L
*
            C = DNRM2( L-K+1, A( K, I ), 1 )
            R = DNRM2( L-K+1, A( I, K ), LDA )
            ICA = IDAMAX( L, A( 1, I ), 1 )
            CA = ABS( A( ICA, I ) )
            IRA = IDAMAX( N-K+1, A( I, K ), LDA )
            RA = ABS( A( I, IRA+K-1 ) )
*
*           Guard against zero C or R due to underflow.
*
            IF( C.EQ.ZERO .OR. R.EQ.ZERO ) CYCLE
*
*           Exit if NaN to avoid infinite loop
*
            IF( DISNAN( C+CA+R+RA ) ) THEN
               INFO = -3
               CALL XERBLA( 'DGEBAL', -INFO )
               RETURN
            END IF
*
            G = R / SCLFAC
            F = ONE
            S = C + R
*
            DO WHILE( C.LT.G .AND. MAX( F, C, CA ).LT.SFMAX2 .AND.
     $                MIN( R, G, RA ).GT.SFMIN2 )
               F = F*SCLFAC
               C = C*SCLFAC
               CA = CA*SCLFAC
               R = R / SCLFAC
               G = G / SCLFAC
               RA = RA / SCLFAC
            END DO
*
            G = C / SCLFAC
*
            DO WHILE( G.GE.R .AND. MAX( R, RA ).LT.SFMAX2 .AND.
     $                MIN( F, C, G, CA ).GT.SFMIN2 )
               F = F / SCLFAC
               C = C / SCLFAC
               G = G / SCLFAC
               CA = CA / SCLFAC
               R = R*SCLFAC
               RA = RA*SCLFAC
            END DO
*
*           Now balance.
*
            IF( ( C+R ).GE.FACTOR*S ) CYCLE
            IF( F.LT.ONE .AND. SCALE( I ).LT.ONE ) THEN
               IF( F*SCALE( I ).LE.SFMIN1 ) CYCLE
            END IF
            IF( F.GT.ONE .AND. SCALE( I ).GT.ONE ) THEN
               IF( SCALE( I ).GE.SFMAX1 / F ) CYCLE
            END IF
            G = ONE / F
            SCALE( I ) = SCALE( I )*F
            NOCONV = .TRUE.
*
            CALL DSCAL( N-K+1, G, A( I, K ), LDA )
            CALL DSCAL( L, F, A( 1, I ), 1 )
*
         END DO
*
      END DO
*
      ILO = K
      IHI = L
*
      RETURN
*
*     End of DGEBAL
*
      END