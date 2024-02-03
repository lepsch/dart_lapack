      SUBROUTINE CGBSVX( FACT, TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, KL, KU, LDAB, LDAFB, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               BERR( * ), C( * ), FERR( * ), R( * ), RWORK( * )       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================
*  Moved setting of INFO = N+1 so INFO does not subsequently get
*  overwritten.  Sven, 17 Mar 05.
*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      String             NORM;
      int                I, INFEQU, J, J1, J2;
      REAL               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               CLANGB, CLANTB, SLAMCH
      // EXTERNAL LSAME, CLANGB, CLANTB, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGBCON, CGBEQU, CGBRFS, CGBTRF, CGBTRS, CLACPY, CLAQGB, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      if ( NOFACT .OR. EQUIL ) {
         EQUED = 'N'
         ROWEQU = .FALSE.
         COLEQU = .FALSE.
      } else {
         ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
         COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         SMLNUM = SLAMCH( 'Safe minimum' )
         BIGNUM = ONE / SMLNUM
      }

      // Test the input parameters.

      if ( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.LSAME( FACT, 'F' ) ) {
         INFO = -1
      } else if ( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( KL.LT.0 ) {
         INFO = -4
      } else if ( KU.LT.0 ) {
         INFO = -5
      } else if ( NRHS.LT.0 ) {
         INFO = -6
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = -8
      } else if ( LDAFB.LT.2*KL+KU+1 ) {
         INFO = -10
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -12
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 10 J = 1, N
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
   10       CONTINUE
            if ( RCMIN.LE.ZERO ) {
               INFO = -13
            } else if ( N.GT.0 ) {
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               ROWCND = ONE
            }
         }
         if ( COLEQU .AND. INFO.EQ.0 ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            DO 20 J = 1, N
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
   20       CONTINUE
            if ( RCMIN.LE.ZERO ) {
               INFO = -14
            } else if ( N.GT.0 ) {
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            }
         }
         if ( INFO.EQ.0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -16
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -18
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('CGBSVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         cgbequ(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU.EQ.0 ) {

            // Equilibrate the matrix.

            claqgb(N, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         }
      }

      // Scale the right hand side.

      if ( NOTRAN ) {
         if ( ROWEQU ) {
            DO 40 J = 1, NRHS
               DO 30 I = 1, N
                  B( I, J ) = R( I )*B( I, J )
   30          CONTINUE
   40       CONTINUE
         }
      } else if ( COLEQU ) {
         DO 60 J = 1, NRHS
            DO 50 I = 1, N
               B( I, J ) = C( I )*B( I, J )
   50       CONTINUE
   60    CONTINUE
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the LU factorization of the band matrix A.

         DO 70 J = 1, N
            J1 = MAX( J-KU, 1 )
            J2 = MIN( J+KL, N )
            ccopy(J2-J1+1, AB( KU+1-J+J1, J ), 1, AFB( KL+KU+1-J+J1, J ), 1 );
   70    CONTINUE

         cgbtrf(N, N, KL, KU, AFB, LDAFB, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {

            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            ANORM = ZERO
            DO 90 J = 1, INFO
               DO 80 I = MAX( KU+2-J, 1 ), MIN( N+KU+1-J, KL+KU+1 )
                  ANORM = MAX( ANORM, ABS( AB( I, J ) ) )
   80          CONTINUE
   90       CONTINUE
            RPVGRW = CLANTB( 'M', 'U', 'N', INFO, MIN( INFO-1, KL+KU ), AFB( MAX( 1, KL+KU+2-INFO ), 1 ), LDAFB, RWORK )
            if ( RPVGRW.EQ.ZERO ) {
               RPVGRW = ONE
            } else {
               RPVGRW = ANORM / RPVGRW
            }
            RWORK( 1 ) = RPVGRW
            RCOND = ZERO
            RETURN
         }
      }

      // Compute the norm of the matrix A and the
      // reciprocal pivot growth factor RPVGRW.

      if ( NOTRAN ) {
         NORM = '1'
      } else {
         NORM = 'I'
      }
      ANORM = CLANGB( NORM, N, KL, KU, AB, LDAB, RWORK )
      RPVGRW = CLANTB( 'M', 'U', 'N', N, KL+KU, AFB, LDAFB, RWORK )
      if ( RPVGRW.EQ.ZERO ) {
         RPVGRW = ONE
      } else {
         RPVGRW = CLANGB( 'M', N, KL, KU, AB, LDAB, RWORK ) / RPVGRW
      }

      // Compute the reciprocal of the condition number of A.

      cgbcon(NORM, N, KL, KU, AFB, LDAFB, IPIV, ANORM, RCOND, WORK, RWORK, INFO );

      // Compute the solution matrix X.

      clacpy('Full', N, NRHS, B, LDB, X, LDX );
      cgbtrs(TRANS, N, KL, KU, NRHS, AFB, LDAFB, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      cgbrfs(TRANS, N, KL, KU, NRHS, AB, LDAB, AFB, LDAFB, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, RWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( NOTRAN ) {
         if ( COLEQU ) {
            DO 110 J = 1, NRHS
               DO 100 I = 1, N
                  X( I, J ) = C( I )*X( I, J )
  100          CONTINUE
  110       CONTINUE
            DO 120 J = 1, NRHS
               FERR( J ) = FERR( J ) / COLCND
  120       CONTINUE
         }
      } else if ( ROWEQU ) {
         DO 140 J = 1, NRHS
            DO 130 I = 1, N
               X( I, J ) = R( I )*X( I, J )
  130       CONTINUE
  140    CONTINUE
         DO 150 J = 1, NRHS
            FERR( J ) = FERR( J ) / ROWCND
  150    CONTINUE
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1

      RWORK( 1 ) = RPVGRW
      RETURN

      // End of CGBSVX

      }
