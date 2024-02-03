      SUBROUTINE ZLATM4( ITYPE, N, NZ1, NZ2, RSIGN, AMAGN, RCOND, TRIANG, IDIST, ISEED, A, LDA )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      bool               RSIGN;
      int                IDIST, ITYPE, LDA, N, NZ1, NZ2;
      double             AMAGN, RCOND, TRIANG;
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      int                I, ISDB, ISDE, JC, JD, JR, K, KBEG, KEND, KLEN;
      double             ALPHA;
      COMPLEX*16         CTEMP
*     ..
*     .. External Functions ..
      double             DLARAN;
      COMPLEX*16         ZLARND
      EXTERNAL           DLARAN, ZLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASET
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, EXP, LOG, MAX, MIN, MOD
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) RETURN
      CALL ZLASET( 'Full', N, N, CZERO, CZERO, A, LDA )
*
*     Insure a correct ISEED
*
      IF( MOD( ISEED( 4 ), 2 ).NE.1 ) ISEED( 4 ) = ISEED( 4 ) + 1
*
*     Compute diagonal and subdiagonal according to ITYPE, NZ1, NZ2,
*     and RCOND
*
      IF( ITYPE.NE.0 ) THEN
         IF( ABS( ITYPE ).GE.4 ) THEN
            KBEG = MAX( 1, MIN( N, NZ1+1 ) )
            KEND = MAX( KBEG, MIN( N, N-NZ2 ) )
            KLEN = KEND + 1 - KBEG
         ELSE
            KBEG = 1
            KEND = N
            KLEN = N
         END IF
         ISDB = 1
         ISDE = 0
         GO TO ( 10, 30, 50, 80, 100, 120, 140, 160, 180, 200 )ABS( ITYPE )
*
*        abs(ITYPE) = 1: Identity
*
   10    CONTINUE
         DO 20 JD = 1, N
            A( JD, JD ) = CONE
   20    CONTINUE
         GO TO 220
*
*        abs(ITYPE) = 2: Transposed Jordan block
*
   30    CONTINUE
         DO 40 JD = 1, N - 1
            A( JD+1, JD ) = CONE
   40    CONTINUE
         ISDB = 1
         ISDE = N - 1
         GO TO 220
*
*        abs(ITYPE) = 3: Transposed Jordan block, followed by the
*                        identity.
*
   50    CONTINUE
         K = ( N-1 ) / 2
         DO 60 JD = 1, K
            A( JD+1, JD ) = CONE
   60    CONTINUE
         ISDB = 1
         ISDE = K
         DO 70 JD = K + 2, 2*K + 1
            A( JD, JD ) = CONE
   70    CONTINUE
         GO TO 220
*
*        abs(ITYPE) = 4: 1,...,k
*
   80    CONTINUE
         DO 90 JD = KBEG, KEND
            A( JD, JD ) = DCMPLX( JD-NZ1 )
   90    CONTINUE
         GO TO 220
*
*        abs(ITYPE) = 5: One large D value:
*
  100    CONTINUE
         DO 110 JD = KBEG + 1, KEND
            A( JD, JD ) = DCMPLX( RCOND )
  110    CONTINUE
         A( KBEG, KBEG ) = CONE
         GO TO 220
*
*        abs(ITYPE) = 6: One small D value:
*
  120    CONTINUE
         DO 130 JD = KBEG, KEND - 1
            A( JD, JD ) = CONE
  130    CONTINUE
         A( KEND, KEND ) = DCMPLX( RCOND )
         GO TO 220
*
*        abs(ITYPE) = 7: Exponentially distributed D values:
*
  140    CONTINUE
         A( KBEG, KBEG ) = CONE
         IF( KLEN.GT.1 ) THEN
            ALPHA = RCOND**( ONE / DBLE( KLEN-1 ) )
            DO 150 I = 2, KLEN
               A( NZ1+I, NZ1+I ) = DCMPLX( ALPHA**DBLE( I-1 ) )
  150       CONTINUE
         END IF
         GO TO 220
*
*        abs(ITYPE) = 8: Arithmetically distributed D values:
*
  160    CONTINUE
         A( KBEG, KBEG ) = CONE
         IF( KLEN.GT.1 ) THEN
            ALPHA = ( ONE-RCOND ) / DBLE( KLEN-1 )
            DO 170 I = 2, KLEN
               A( NZ1+I, NZ1+I ) = DCMPLX( DBLE( KLEN-I )*ALPHA+RCOND )
  170       CONTINUE
         END IF
         GO TO 220
*
*        abs(ITYPE) = 9: Randomly distributed D values on ( RCOND, 1):
*
  180    CONTINUE
         ALPHA = LOG( RCOND )
         DO 190 JD = KBEG, KEND
            A( JD, JD ) = EXP( ALPHA*DLARAN( ISEED ) )
  190    CONTINUE
         GO TO 220
*
*        abs(ITYPE) = 10: Randomly distributed D values from DIST
*
  200    CONTINUE
         DO 210 JD = KBEG, KEND
            A( JD, JD ) = ZLARND( IDIST, ISEED )
  210    CONTINUE
*
  220    CONTINUE
*
*        Scale by AMAGN
*
         DO 230 JD = KBEG, KEND
            A( JD, JD ) = AMAGN*DBLE( A( JD, JD ) )
  230    CONTINUE
         DO 240 JD = ISDB, ISDE
            A( JD+1, JD ) = AMAGN*DBLE( A( JD+1, JD ) )
  240    CONTINUE
*
*        If RSIGN = .TRUE., assign random signs to diagonal and
*        subdiagonal
*
         IF( RSIGN ) THEN
            DO 250 JD = KBEG, KEND
               IF( DBLE( A( JD, JD ) ).NE.ZERO ) THEN
                  CTEMP = ZLARND( 3, ISEED )
                  CTEMP = CTEMP / ABS( CTEMP )
                  A( JD, JD ) = CTEMP*DBLE( A( JD, JD ) )
               END IF
  250       CONTINUE
            DO 260 JD = ISDB, ISDE
               IF( DBLE( A( JD+1, JD ) ).NE.ZERO ) THEN
                  CTEMP = ZLARND( 3, ISEED )
                  CTEMP = CTEMP / ABS( CTEMP )
                  A( JD+1, JD ) = CTEMP*DBLE( A( JD+1, JD ) )
               END IF
  260       CONTINUE
         END IF
*
*        Reverse if ITYPE < 0
*
         IF( ITYPE.LT.0 ) THEN
            DO 270 JD = KBEG, ( KBEG+KEND-1 ) / 2
               CTEMP = A( JD, JD )
               A( JD, JD ) = A( KBEG+KEND-JD, KBEG+KEND-JD )
               A( KBEG+KEND-JD, KBEG+KEND-JD ) = CTEMP
  270       CONTINUE
            DO 280 JD = 1, ( N-1 ) / 2
               CTEMP = A( JD+1, JD )
               A( JD+1, JD ) = A( N+1-JD, N-JD )
               A( N+1-JD, N-JD ) = CTEMP
  280       CONTINUE
         END IF
*
      END IF
*
*     Fill in upper triangle
*
      IF( TRIANG.NE.ZERO ) THEN
         DO 300 JC = 2, N
            DO 290 JR = 1, JC - 1
               A( JR, JC ) = TRIANG*ZLARND( IDIST, ISEED )
  290       CONTINUE
  300    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLATM4
*
      END
