      SUBROUTINE ZGSVJ1( JOBV, M, N, N1, A, LDA, D, SVA, MV, V, LDV, EPS, SFMIN, TOL, NSWEEP, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*     .. Scalar Arguments ..
      DOUBLE PRECISION   EPS, SFMIN, TOL
      int                INFO, LDA, LDV, LWORK, M, MV, N, N1, NSWEEP
      String             JOBV;
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), D( N ), V( LDV, * ), WORK( LWORK )
      DOUBLE PRECISION   SVA( N )
*     ..
*
*  =====================================================================
*
*     .. Local Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*     ..
*     .. Local Scalars ..
      COMPLEX*16         AAPQ, OMPQ
      DOUBLE PRECISION   AAPP, AAPP0, AAPQ1, AAQQ, APOAQ, AQOAP, BIG, BIGTHETA, CS, MXAAPQ, MXSINJ, ROOTBIG, ROOTEPS, ROOTSFMIN, ROOTTOL, SMALL, SN, T, TEMP1, THETA, THSIGN
      int                BLSKIP, EMPTSW, i, ibr, igl, IERR, IJBLSK, ISWROT, jbc, jgl, KBL, MVL, NOTROT, nblc, nblr, p, PSKIPPED, q, ROWSKIP, SWBAND
      LOGICAL            APPLV, ROTOK, RSVEC
*     ..
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, MAX, DBLE, MIN, SIGN, SQRT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DZNRM2
      COMPLEX*16         ZDOTC
      int                IDAMAX
      LOGICAL            LSAME
      EXTERNAL           IDAMAX, LSAME, ZDOTC, DZNRM2
*     ..
*     .. External Subroutines ..
*     .. from BLAS
      EXTERNAL           ZCOPY, ZROT, ZSWAP, ZAXPY
*     .. from LAPACK
      EXTERNAL           ZLASCL, ZLASSQ, XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      APPLV = LSAME( JOBV, 'A' )
      RSVEC = LSAME( JOBV, 'V' )
      IF( .NOT.( RSVEC .OR. APPLV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( N.LT.0 ) .OR. ( N.GT.M ) ) THEN
         INFO = -3
      ELSE IF( N1.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.M ) THEN
         INFO = -6
      ELSE IF( ( RSVEC.OR.APPLV ) .AND. ( MV.LT.0 ) ) THEN
         INFO = -9
      ELSE IF( ( RSVEC.AND.( LDV.LT.N ) ).OR. ( APPLV.AND.( LDV.LT.MV ) )  ) THEN
         INFO = -11
      ELSE IF( TOL.LE.EPS ) THEN
         INFO = -14
      ELSE IF( NSWEEP.LT.0 ) THEN
         INFO = -15
      ELSE IF( LWORK.LT.M ) THEN
         INFO = -17
      ELSE
         INFO = 0
      END IF
*
*     #:(
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGSVJ1', -INFO )
         RETURN
      END IF
*
      IF( RSVEC ) THEN
         MVL = N
      ELSE IF( APPLV ) THEN
         MVL = MV
      END IF
      RSVEC = RSVEC .OR. APPLV

      ROOTEPS = SQRT( EPS )
      ROOTSFMIN = SQRT( SFMIN )
      SMALL = SFMIN / EPS
      BIG = ONE / SFMIN
      ROOTBIG = ONE / ROOTSFMIN
*     LARGE = BIG / SQRT( DBLE( M*N ) )
      BIGTHETA = ONE / ROOTEPS
      ROOTTOL = SQRT( TOL )
*
*     .. Initialize the right singular vector matrix ..
*
*     RSVEC = LSAME( JOBV, 'Y' )
*
      EMPTSW = N1*( N-N1 )
      NOTROT = 0
*
*     .. Row-cyclic pivot strategy with de Rijk's pivoting ..
*
      KBL = MIN( 8, N )
      NBLR = N1 / KBL
      IF( ( NBLR*KBL ).NE.N1 )NBLR = NBLR + 1

*     .. the tiling is nblr-by-nblc [tiles]

      NBLC = ( N-N1 ) / KBL
      IF( ( NBLC*KBL ).NE.( N-N1 ) )NBLC = NBLC + 1
      BLSKIP = ( KBL**2 ) + 1
*[TP] BLKSKIP is a tuning parameter that depends on SWBAND and KBL.

      ROWSKIP = MIN( 5, KBL )
*[TP] ROWSKIP is a tuning parameter.
      SWBAND = 0
*[TP] SWBAND is a tuning parameter. It is meaningful and effective
*     if ZGESVJ is used as a computational routine in the preconditioned
*     Jacobi SVD algorithm ZGEJSV.
*
*
*     | *   *   * [x] [x] [x]|
*     | *   *   * [x] [x] [x]|    Row-cycling in the nblr-by-nblc [x] blocks.
*     | *   *   * [x] [x] [x]|    Row-cyclic pivoting inside each [x] block.
*     |[x] [x] [x] *   *   * |
*     |[x] [x] [x] *   *   * |
*     |[x] [x] [x] *   *   * |
*
*
      DO 1993 i = 1, NSWEEP
*
*     .. go go go ...
*
         MXAAPQ = ZERO
         MXSINJ = ZERO
         ISWROT = 0
*
         NOTROT = 0
         PSKIPPED = 0
*
*     Each sweep is unrolled using KBL-by-KBL tiles over the pivot pairs
*     1 <= p < q <= N. This is the first step toward a blocked implementation
*     of the rotations. New implementation, based on block transformations,
*     is under development.
*
         DO 2000 ibr = 1, NBLR
*
            igl = ( ibr-1 )*KBL + 1
*

*
* ... go to the off diagonal blocks
*
            igl = ( ibr-1 )*KBL + 1
*
*            DO 2010 jbc = ibr + 1, NBL
            DO 2010 jbc = 1, NBLC
*
               jgl = ( jbc-1 )*KBL + N1 + 1
*
*        doing the block at ( ibr, jbc )
*
               IJBLSK = 0
               DO 2100 p = igl, MIN( igl+KBL-1, N1 )
*
                  AAPP = SVA( p )
                  IF( AAPP.GT.ZERO ) THEN
*
                     PSKIPPED = 0
*
                     DO 2200 q = jgl, MIN( jgl+KBL-1, N )
*
                        AAQQ = SVA( q )
                        IF( AAQQ.GT.ZERO ) THEN
                           AAPP0 = AAPP
*
*     .. M x 2 Jacobi SVD ..
*
*        Safe Gram matrix computation
*
                           IF( AAQQ.GE.ONE ) THEN
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = ( SMALL*AAPP ).LE.AAQQ
                              ELSE
                                 ROTOK = ( SMALL*AAQQ ).LE.AAPP
                              END IF
                              IF( AAPP.LT.( BIG / AAQQ ) ) THEN
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / AAQQ ) / AAPP
                              ELSE
                                 CALL ZCOPY( M, A( 1, p ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK, LDA, IERR )
                                 AAPQ = ZDOTC( M, WORK, 1, A( 1, q ), 1 ) / AAQQ
                              END IF
                           ELSE
                              IF( AAPP.GE.AAQQ ) THEN
                                 ROTOK = AAPP.LE.( AAQQ / SMALL )
                              ELSE
                                 ROTOK = AAQQ.LE.( AAPP / SMALL )
                              END IF
                              IF( AAPP.GT.( SMALL / AAQQ ) ) THEN
                                 AAPQ = ( ZDOTC( M, A( 1, p ), 1, A( 1, q ), 1 ) / MAX(AAQQ,AAPP) ) / MIN(AAQQ,AAPP)
                              ELSE
                                 CALL ZCOPY( M, A( 1, q ), 1, WORK, 1 )                                  CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK, LDA, IERR )
                                 AAPQ = ZDOTC( M, A( 1, p ), 1, WORK, 1 ) / AAPP
                              END IF
                           END IF
*
*                           AAPQ = AAPQ * CONJG(CWORK(p))*CWORK(q)
                           AAPQ1  = -ABS(AAPQ)
                           MXAAPQ = MAX( MXAAPQ, -AAPQ1 )
*
*        TO rotate or NOT to rotate, THAT is the question ...
*
                           IF( ABS( AAPQ1 ).GT.TOL ) THEN
                              OMPQ = AAPQ / ABS(AAPQ)
                              NOTROT = 0
*[RTD]      ROTATED  = ROTATED + 1
                              PSKIPPED = 0
                              ISWROT = ISWROT + 1
*
                              IF( ROTOK ) THEN
*
                                 AQOAP = AAQQ / AAPP
                                 APOAQ = AAPP / AAQQ
                                 THETA = -HALF*ABS( AQOAP-APOAQ )/ AAPQ1
                                 IF( AAQQ.GT.AAPP0 )THETA = -THETA
*
                                 IF( ABS( THETA ).GT.BIGTHETA ) THEN
                                    T  = HALF / THETA
                                    CS = ONE
                                    CALL ZROT( M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*T )
                                    IF( RSVEC ) THEN
                                        CALL ZROT( MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*T )
                                    END IF
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, ABS( T ) )
                                 ELSE
*
*                 .. choose correct signum for THETA and rotate
*
                                    THSIGN = -SIGN( ONE, AAPQ1 )
                                    IF( AAQQ.GT.AAPP0 )THSIGN = -THSIGN
                                    T = ONE / ( THETA+THSIGN* SQRT( ONE+THETA*THETA ) )
                                    CS = SQRT( ONE / ( ONE+T*T ) )
                                    SN = T*CS
                                    MXSINJ = MAX( MXSINJ, ABS( SN ) )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE+T*APOAQ*AAPQ1 ) )                                     AAPP = AAPP*SQRT( MAX( ZERO, ONE-T*AQOAP*AAPQ1 ) )
*
                                    CALL ZROT( M, A(1,p), 1, A(1,q), 1, CS, CONJG(OMPQ)*SN )
                                    IF( RSVEC ) THEN
                                        CALL ZROT( MVL, V(1,p), 1, V(1,q), 1, CS, CONJG(OMPQ)*SN )
                                    END IF
                                 END IF
                                 D(p) = -D(q) * OMPQ
*
                              ELSE
*              .. have to use modified Gram-Schmidt like transformation
                               IF( AAPP.GT.AAQQ ) THEN
                                    CALL ZCOPY( M, A( 1, p ), 1, WORK, 1 )                                     CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, WORK,LDA, IERR )                                     CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, A( 1, q ), LDA, IERR )                                     CALL ZAXPY( M, -AAPQ, WORK, 1, A( 1, q ), 1 )                                     CALL ZLASCL( 'G', 0, 0, ONE, AAQQ, M, 1, A( 1, q ), LDA, IERR )
                                    SVA( q ) = AAQQ*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               ELSE
                                   CALL ZCOPY( M, A( 1, q ), 1, WORK, 1 )                                     CALL ZLASCL( 'G', 0, 0, AAQQ, ONE, M, 1, WORK,LDA, IERR )                                     CALL ZLASCL( 'G', 0, 0, AAPP, ONE, M, 1, A( 1, p ), LDA, IERR )                                     CALL ZAXPY( M, -CONJG(AAPQ), WORK, 1, A( 1, p ), 1 )                                     CALL ZLASCL( 'G', 0, 0, ONE, AAPP, M, 1, A( 1, p ), LDA, IERR )
                                    SVA( p ) = AAPP*SQRT( MAX( ZERO, ONE-AAPQ1*AAPQ1 ) )
                                    MXSINJ = MAX( MXSINJ, SFMIN )
                               END IF
                              END IF
*           END IF ROTOK THEN ... ELSE
*
*           In the case of cancellation in updating SVA(q), SVA(p)
*           .. recompute SVA(q), SVA(p)
                              IF( ( SVA( q ) / AAQQ )**2.LE.ROOTEPS ) THEN                                  IF( ( AAQQ.LT.ROOTBIG ) .AND. ( AAQQ.GT.ROOTSFMIN ) ) THEN
                                    SVA( q ) = DZNRM2( M, A( 1, q ), 1)
                                  ELSE
                                    T = ZERO
                                    AAQQ = ONE
                                    CALL ZLASSQ( M, A( 1, q ), 1, T, AAQQ )
                                    SVA( q ) = T*SQRT( AAQQ )
                                 END IF
                              END IF
                              IF( ( AAPP / AAPP0 )**2.LE.ROOTEPS ) THEN
                                 IF( ( AAPP.LT.ROOTBIG ) .AND. ( AAPP.GT.ROOTSFMIN ) ) THEN
                                    AAPP = DZNRM2( M, A( 1, p ), 1 )
                                 ELSE
                                    T = ZERO
                                    AAPP = ONE
                                    CALL ZLASSQ( M, A( 1, p ), 1, T, AAPP )
                                    AAPP = T*SQRT( AAPP )
                                 END IF
                                 SVA( p ) = AAPP
                              END IF
*              end of OK rotation
                           ELSE
                              NOTROT = NOTROT + 1
*[RTD]      SKIPPED  = SKIPPED  + 1
                              PSKIPPED = PSKIPPED + 1
                              IJBLSK = IJBLSK + 1
                           END IF
                        ELSE
                           NOTROT = NOTROT + 1
                           PSKIPPED = PSKIPPED + 1
                           IJBLSK = IJBLSK + 1
                        END IF
*
                        IF( ( i.LE.SWBAND ) .AND. ( IJBLSK.GE.BLSKIP ) ) THEN
                           SVA( p ) = AAPP
                           NOTROT = 0
                           GO TO 2011
                        END IF
                        IF( ( i.LE.SWBAND ) .AND. ( PSKIPPED.GT.ROWSKIP ) ) THEN
                           AAPP = -AAPP
                           NOTROT = 0
                           GO TO 2203
                        END IF
*
 2200                CONTINUE
*        end of the q-loop
 2203                CONTINUE
*
                     SVA( p ) = AAPP
*
                  ELSE
*
                     IF( AAPP.EQ.ZERO )NOTROT = NOTROT + MIN( jgl+KBL-1, N ) - jgl + 1
                     IF( AAPP.LT.ZERO )NOTROT = 0
*
                  END IF
*
 2100          CONTINUE
*     end of the p-loop
 2010       CONTINUE
*     end of the jbc-loop
 2011       CONTINUE
*2011 bailed out of the jbc-loop
            DO 2012 p = igl, MIN( igl+KBL-1, N )
               SVA( p ) = ABS( SVA( p ) )
 2012       CONTINUE
***
 2000    CONTINUE
*2000 :: end of the ibr-loop
*
*     .. update SVA(N)
         IF( ( SVA( N ).LT.ROOTBIG ) .AND. ( SVA( N ).GT.ROOTSFMIN ) ) THEN
            SVA( N ) = DZNRM2( M, A( 1, N ), 1 )
         ELSE
            T = ZERO
            AAPP = ONE
            CALL ZLASSQ( M, A( 1, N ), 1, T, AAPP )
            SVA( N ) = T*SQRT( AAPP )
         END IF
*
*     Additional steering devices
*
         IF( ( i.LT.SWBAND ) .AND. ( ( MXAAPQ.LE.ROOTTOL ) .OR. ( ISWROT.LE.N ) ) )SWBAND = i
*
         IF( ( i.GT.SWBAND+1 ) .AND. ( MXAAPQ.LT.SQRT( DBLE( N ) )* TOL ) .AND. ( DBLE( N )*MXAAPQ*MXSINJ.LT.TOL ) ) THEN
            GO TO 1994
         END IF
*
         IF( NOTROT.GE.EMPTSW )GO TO 1994
*
 1993 CONTINUE
*     end i=1:NSWEEP loop
*
* #:( Reaching this point means that the procedure has not converged.
      INFO = NSWEEP - 1
      GO TO 1995
*
 1994 CONTINUE
* #:) Reaching this point means numerical convergence after the i-th
*     sweep.
*
      INFO = 0
* #:) INFO = 0 confirms successful iterations.
 1995 CONTINUE
*
*     Sort the vector SVA() of column norms.
      DO 5991 p = 1, N - 1
         q = IDAMAX( N-p+1, SVA( p ), 1 ) + p - 1
         IF( p.NE.q ) THEN
            TEMP1 = SVA( p )
            SVA( p ) = SVA( q )
            SVA( q ) = TEMP1
            AAPQ = D( p )
            D( p ) = D( q )
            D( q ) = AAPQ
            CALL ZSWAP( M, A( 1, p ), 1, A( 1, q ), 1 )
            IF( RSVEC )CALL ZSWAP( MVL, V( 1, p ), 1, V( 1, q ), 1 )
         END IF
 5991 CONTINUE
*
*
      RETURN
*     ..
*     .. END OF ZGSVJ1
*     ..
      END
