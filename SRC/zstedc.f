      SUBROUTINE ZSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, LIWORK, LRWORK, LWORK, N;
*     ..
*     .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), RWORK( * );
      COMPLEX*16         WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE, TWO;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      bool               LQUERY;
      int                FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL, LRWMIN, LWMIN, M, SMLSIZ, START;
      double             EPS, ORGNRM, P, TINY;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANST
*     ..
*     .. External Subroutines ..
      // EXTERNAL DLASCL, DLASET, DSTEDC, DSTEQR, DSTERF, XERBLA, ZLACPY, ZLACRM, ZLAED0, ZSTEQR, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX, MOD, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LRWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Compute the workspace requirements
*
         SMLSIZ = ILAENV( 9, 'ZSTEDC', ' ', 0, 0, 0, 0 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LWMIN = 1
            LIWMIN = 1
            LRWMIN = 2*( N - 1 )
         ELSE IF( ICOMPZ.EQ.1 ) THEN
            LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1
            LWMIN = N*N
            LRWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWMIN = 6 + 6*N + 5*N*LGN
         ELSE IF( ICOMPZ.EQ.2 ) THEN
            LWMIN = 1
            LRWMIN = 1 + 4*N + 2*N**2
            LIWMIN = 3 + 5*N
         END IF
         WORK( 1 ) = LWMIN
         RWORK( 1 ) = LRWMIN
         IWORK( 1 ) = LIWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -8
         ELSE IF( LRWORK.LT.LRWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -10
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEDC', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     If the following conditional clause is removed, then the routine
*     will use the Divide and Conquer routine to compute only the
*     eigenvalues, which requires (3N + 3N**2) real workspace and
*     (2 + 5N + 2N lg(N)) integer workspace.
*     Since on many architectures DSTERF is much faster than any other
*     algorithm for finding eigenvalues only, it is used here
*     as the default. If the conditional clause is removed, then
*     information on the size of workspace needs to be changed.
*
*     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
*
      IF( ICOMPZ.EQ.0 ) THEN
         CALL DSTERF( N, D, E, INFO )
         GO TO 70
      END IF
*
*     If N is smaller than the minimum divide size (SMLSIZ+1), then
*     solve the problem with another solver.
*
      IF( N.LE.SMLSIZ ) THEN
*
         CALL ZSTEQR( COMPZ, N, D, E, Z, LDZ, RWORK, INFO )
*
      ELSE
*
*        If COMPZ = 'I', we simply call DSTEDC instead.
*
         IF( ICOMPZ.EQ.2 ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, RWORK, N )
            LL = N*N + 1
            CALL DSTEDC( 'I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1, IWORK, LIWORK, INFO )
            DO 20 J = 1, N
               DO 10 I = 1, N
                  Z( I, J ) = RWORK( ( J-1 )*N+I )
   10          CONTINUE
   20       CONTINUE
            GO TO 70
         END IF
*
*        From now on, only option left to be handled is COMPZ = 'V',
*        i.e. ICOMPZ = 1.
*
*        Scale.
*
         ORGNRM = DLANST( 'M', N, D, E )
         IF( ORGNRM.EQ.ZERO ) GO TO 70
*
         EPS = DLAMCH( 'Epsilon' )
*
         START = 1
*
*        while ( START <= N )
*
   30    CONTINUE
         IF( START.LE.N ) THEN
*
*           Let FINISH be the position of the next subdiagonal entry
*           such that E( FINISH ) <= TINY or FINISH = N if no such
*           subdiagonal exists.  The matrix identified by the elements
*           between START and FINISH constitutes an independent
*           sub-problem.
*
            FINISH = START
   40       CONTINUE
            IF( FINISH.LT.N ) THEN
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 40
               END IF
            END IF
*
*           (Sub) Problem determined.  Compute its size and solve it.
*
            M = FINISH - START + 1
            IF( M.GT.SMLSIZ ) THEN
*
*              Scale.
*
               ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
               CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO )                CALL DLASCL( 'G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO )
*
               CALL ZLAED0( N, M, D( START ), E( START ), Z( 1, START ), LDZ, WORK, N, RWORK, IWORK, INFO )
               IF( INFO.GT.0 ) THEN
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 70
               END IF
*
*              Scale back.
*
               CALL DLASCL( 'G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO )
*
            ELSE
               CALL DSTEQR( 'I', M, D( START ), E( START ), RWORK, M, RWORK( M*M+1 ), INFO )                CALL ZLACRM( N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N, RWORK( M*M+1 ) )
               CALL ZLACPY( 'A', N, M, WORK, N, Z( 1, START ), LDZ )
               IF( INFO.GT.0 ) THEN
                  INFO = START*( N+1 ) + FINISH
                  GO TO 70
               END IF
            END IF
*
            START = FINISH + 1
            GO TO 30
         END IF
*
*        endwhile
*
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 60 II = 2, N
           I = II - 1
           K = I
           P = D( I )
           DO 50 J = II, N
              IF( D( J ).LT.P ) THEN
                 K = J
                 P = D( J )
              END IF
   50      CONTINUE
           IF( K.NE.I ) THEN
              D( K ) = D( I )
              D( I ) = P
              CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
           END IF
   60    CONTINUE
      END IF
*
   70 CONTINUE
      WORK( 1 ) = LWMIN
      RWORK( 1 ) = LRWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of ZSTEDC
*
      END
