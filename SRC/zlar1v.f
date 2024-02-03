      SUBROUTINE ZLAR1V( N, B1, BN, LAMBDA, D, L, LD, LLD, PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ, MINGMA, R, ISUPPZ, NRMINV, RESID, RQCORR, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      bool               WANTNC;
      int       B1, BN, N, NEGCNT, R;
      double             GAPTOL, LAMBDA, MINGMA, NRMINV, PIVMIN, RESID, RQCORR, ZTZ;
*     ..
*     .. Array Arguments ..
      int                ISUPPZ( * );
      double             D( * ), L( * ), LD( * ), LLD( * ), WORK( * );
      COMPLEX*16       Z( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )

*     ..
*     .. Local Scalars ..
      bool               SAWNAN1, SAWNAN2;
      int                I, INDLPL, INDP, INDS, INDUMN, NEG1, NEG2, R1, R2;
      double             DMINUS, DPLUS, EPS, S, TMP;
*     ..
*     .. External Functions ..
      bool    DISNAN;
      double             DLAMCH;
      EXTERNAL           DISNAN, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE
*     ..
*     .. Executable Statements ..
*
      EPS = DLAMCH( 'Precision' )


      IF( R.EQ.0 ) THEN
         R1 = B1
         R2 = BN
      ELSE
         R1 = R
         R2 = R
      END IF

*     Storage for LPLUS
      INDLPL = 0
*     Storage for UMINUS
      INDUMN = N
      INDS = 2*N + 1
      INDP = 3*N + 1

      IF( B1.EQ.1 ) THEN
         WORK( INDS ) = ZERO
      ELSE
         WORK( INDS+B1-1 ) = LLD( B1-1 )
      END IF

*
*     Compute the stationary transform (using the differential form)
*     until the index R2.
*
      SAWNAN1 = .FALSE.
      NEG1 = 0
      S = WORK( INDS+B1-1 ) - LAMBDA
      DO 50 I = B1, R1 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 50   CONTINUE
      SAWNAN1 = DISNAN( S )
      IF( SAWNAN1 ) GOTO 60
      DO 51 I = R1, R2 - 1
         DPLUS = D( I ) + S
         WORK( INDLPL+I ) = LD( I ) / DPLUS
         WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
         S = WORK( INDS+I ) - LAMBDA
 51   CONTINUE
      SAWNAN1 = DISNAN( S )
*
 60   CONTINUE
      IF( SAWNAN1 ) THEN
*        Runs a slower version of the above loop if a NaN is detected
         NEG1 = 0
         S = WORK( INDS+B1-1 ) - LAMBDA
         DO 70 I = B1, R1 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            IF(DPLUS.LT.ZERO) NEG1 = NEG1 + 1
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO ) WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 70      CONTINUE
         DO 71 I = R1, R2 - 1
            DPLUS = D( I ) + S
            IF(ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK( INDLPL+I ) = LD( I ) / DPLUS
            WORK( INDS+I ) = S*WORK( INDLPL+I )*L( I )
            IF( WORK( INDLPL+I ).EQ.ZERO ) WORK( INDS+I ) = LLD( I )
            S = WORK( INDS+I ) - LAMBDA
 71      CONTINUE
      END IF
*
*     Compute the progressive transform (using the differential form)
*     until the index R1
*
      SAWNAN2 = .FALSE.
      NEG2 = 0
      WORK( INDP+BN-1 ) = D( BN ) - LAMBDA
      DO 80 I = BN - 1, R1, -1
         DMINUS = LLD( I ) + WORK( INDP+I )
         TMP = D( I ) / DMINUS
         IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
         WORK( INDUMN+I ) = L( I )*TMP
         WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
 80   CONTINUE
      TMP = WORK( INDP+R1-1 )
      SAWNAN2 = DISNAN( TMP )

      IF( SAWNAN2 ) THEN
*        Runs a slower version of the above loop if a NaN is detected
         NEG2 = 0
         DO 100 I = BN-1, R1, -1
            DMINUS = LLD( I ) + WORK( INDP+I )
            IF(ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D( I ) / DMINUS
            IF(DMINUS.LT.ZERO) NEG2 = NEG2 + 1
            WORK( INDUMN+I ) = L( I )*TMP
            WORK( INDP+I-1 ) = WORK( INDP+I )*TMP - LAMBDA
            IF( TMP.EQ.ZERO ) WORK( INDP+I-1 ) = D( I ) - LAMBDA
 100     CONTINUE
      END IF
*
*     Find the index (from R1 to R2) of the largest (in magnitude)
*     diagonal element of the inverse
*
      MINGMA = WORK( INDS+R1-1 ) + WORK( INDP+R1-1 )
      IF( MINGMA.LT.ZERO ) NEG1 = NEG1 + 1
      IF( WANTNC ) THEN
         NEGCNT = NEG1 + NEG2
      ELSE
         NEGCNT = -1
      ENDIF
      IF( ABS(MINGMA).EQ.ZERO ) MINGMA = EPS*WORK( INDS+R1-1 )
      R = R1
      DO 110 I = R1, R2 - 1
         TMP = WORK( INDS+I ) + WORK( INDP+I )
         IF( TMP.EQ.ZERO ) TMP = EPS*WORK( INDS+I )
         IF( ABS( TMP ).LE.ABS( MINGMA ) ) THEN
            MINGMA = TMP
            R = I + 1
         END IF
 110  CONTINUE
*
*     Compute the FP vector: solve N^T v = e_r
*
      ISUPPZ( 1 ) = B1
      ISUPPZ( 2 ) = BN
      Z( R ) = CONE
      ZTZ = ONE
*
*     Compute the FP vector upwards from R
*
      IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
         DO 210 I = R-1, B1, -1
            Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL ) THEN
               Z( I ) = ZERO
               ISUPPZ( 1 ) = I + 1
               GOTO 220
            ENDIF
            ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
 210     CONTINUE
 220     CONTINUE
      ELSE
*        Run slower loop if NaN occurred.
         DO 230 I = R - 1, B1, -1
            IF( Z( I+1 ).EQ.ZERO ) THEN
               Z( I ) = -( LD( I+1 ) / LD( I ) )*Z( I+2 )
            ELSE
               Z( I ) = -( WORK( INDLPL+I )*Z( I+1 ) )
            END IF
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL ) THEN
               Z( I ) = ZERO
               ISUPPZ( 1 ) = I + 1
               GO TO 240
            END IF
            ZTZ = ZTZ + DBLE( Z( I )*Z( I ) )
 230     CONTINUE
 240     CONTINUE
      ENDIF

*     Compute the FP vector downwards from R in blocks of size BLKSIZ
      IF( .NOT.SAWNAN1 .AND. .NOT.SAWNAN2 ) THEN
         DO 250 I = R, BN-1
            Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL ) THEN
               Z( I+1 ) = ZERO
               ISUPPZ( 2 ) = I
               GO TO 260
            END IF
            ZTZ = ZTZ + DBLE( Z( I+1 )*Z( I+1 ) )
 250     CONTINUE
 260     CONTINUE
      ELSE
*        Run slower loop if NaN occurred.
         DO 270 I = R, BN - 1
            IF( Z( I ).EQ.ZERO ) THEN
               Z( I+1 ) = -( LD( I-1 ) / LD( I ) )*Z( I-1 )
            ELSE
               Z( I+1 ) = -( WORK( INDUMN+I )*Z( I ) )
            END IF
            IF( (ABS(Z(I))+ABS(Z(I+1)))* ABS(LD(I)).LT.GAPTOL ) THEN
               Z( I+1 ) = ZERO
               ISUPPZ( 2 ) = I
               GO TO 280
            END IF
            ZTZ = ZTZ + DBLE( Z( I+1 )*Z( I+1 ) )
 270     CONTINUE
 280     CONTINUE
      END IF
*
*     Compute quantities for convergence test
*
      TMP = ONE / ZTZ
      NRMINV = SQRT( TMP )
      RESID = ABS( MINGMA )*NRMINV
      RQCORR = MINGMA*TMP
*
*
      RETURN
*
*     End of ZLAR1V
*
      END
