      SUBROUTINE SGESVX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, EQUED, R, C, B, LDB, X, LDX, RCOND, FERR, BERR, WORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, FACT, TRANS;
      int                INFO, LDA, LDAF, LDB, LDX, N, NRHS;
      REAL               RCOND
      // ..
      // .. Array Arguments ..
      int                IPIV( * ), IWORK( * );
      REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ), BERR( * ), C( * ), FERR( * ), R( * ), WORK( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               COLEQU, EQUIL, NOFACT, NOTRAN, ROWEQU;
      String             NORM;
      int                I, INFEQU, J;
      REAL               AMAX, ANORM, BIGNUM, COLCND, RCMAX, RCMIN, ROWCND, RPVGRW, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLANGE, SLANTR
      // EXTERNAL LSAME, SLAMCH, SLANGE, SLANTR
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGECON, SGEEQU, SGERFS, SGETRF, SGETRS, SLACPY, SLAQGE, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      NOTRAN = LSAME( TRANS, 'N' )
      if ( NOFACT .OR. EQUIL ) {
         EQUED = 'N'
         ROWEQU = false;
         COLEQU = false;
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
      } else if ( NRHS.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDAF.LT.MAX( 1, N ) ) {
         INFO = -8
      } else if ( LSAME( FACT, 'F' ) .AND. .NOT. ( ROWEQU .OR. COLEQU .OR. LSAME( EQUED, 'N' ) ) ) {
         INFO = -10
      } else {
         if ( ROWEQU ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 10
               RCMIN = MIN( RCMIN, R( J ) )
               RCMAX = MAX( RCMAX, R( J ) )
            } // 10
            if ( RCMIN.LE.ZERO ) {
               INFO = -11
            } else if ( N.GT.0 ) {
               ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               ROWCND = ONE
            }
         }
         if ( COLEQU .AND. INFO == 0 ) {
            RCMIN = BIGNUM
            RCMAX = ZERO
            for (J = 1; J <= N; J++) { // 20
               RCMIN = MIN( RCMIN, C( J ) )
               RCMAX = MAX( RCMAX, C( J ) )
            } // 20
            if ( RCMIN.LE.ZERO ) {
               INFO = -12
            } else if ( N.GT.0 ) {
               COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
            } else {
               COLCND = ONE
            }
         }
         if ( INFO == 0 ) {
            if ( LDB.LT.MAX( 1, N ) ) {
               INFO = -14
            } else if ( LDX.LT.MAX( 1, N ) ) {
               INFO = -16
            }
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SGESVX', -INFO );
         RETURN
      }

      if ( EQUIL ) {

         // Compute row and column scalings to equilibrate the matrix A.

         sgeequ(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFEQU );
         if ( INFEQU == 0 ) {

            // Equilibrate the matrix.

            slaqge(N, N, A, LDA, R, C, ROWCND, COLCND, AMAX, EQUED );
            ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
            COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
         }
      }

      // Scale the right hand side.

      if ( NOTRAN ) {
         if ( ROWEQU ) {
            for (J = 1; J <= NRHS; J++) { // 40
               for (I = 1; I <= N; I++) { // 30
                  B( I, J ) = R( I )*B( I, J )
               } // 30
            } // 40
         }
      } else if ( COLEQU ) {
         for (J = 1; J <= NRHS; J++) { // 60
            for (I = 1; I <= N; I++) { // 50
               B( I, J ) = C( I )*B( I, J )
            } // 50
         } // 60
      }

      if ( NOFACT .OR. EQUIL ) {

         // Compute the LU factorization of A.

         slacpy('Full', N, N, A, LDA, AF, LDAF );
         sgetrf(N, N, AF, LDAF, IPIV, INFO );

         // Return if INFO is non-zero.

         if ( INFO.GT.0 ) {

            // Compute the reciprocal pivot growth factor of the
            // leading rank-deficient INFO columns of A.

            RPVGRW = SLANTR( 'M', 'U', 'N', INFO, INFO, AF, LDAF, WORK )
            if ( RPVGRW == ZERO ) {
               RPVGRW = ONE
            } else {
               RPVGRW = SLANGE( 'M', N, INFO, A, LDA, WORK ) / RPVGRW
            }
            WORK( 1 ) = RPVGRW
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
      ANORM = SLANGE( NORM, N, N, A, LDA, WORK )
      RPVGRW = SLANTR( 'M', 'U', 'N', N, N, AF, LDAF, WORK )
      if ( RPVGRW == ZERO ) {
         RPVGRW = ONE
      } else {
         RPVGRW = SLANGE( 'M', N, N, A, LDA, WORK ) / RPVGRW
      }

      // Compute the reciprocal of the condition number of A.

      sgecon(NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO );

      // Compute the solution matrix X.

      slacpy('Full', N, NRHS, B, LDB, X, LDX );
      sgetrs(TRANS, N, NRHS, AF, LDAF, IPIV, X, LDX, INFO );

      // Use iterative refinement to improve the computed solution and
      // compute error bounds and backward error estimates for it.

      sgerfs(TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV, B, LDB, X, LDX, FERR, BERR, WORK, IWORK, INFO );

      // Transform the solution matrix X to a solution of the original
      // system.

      if ( NOTRAN ) {
         if ( COLEQU ) {
            for (J = 1; J <= NRHS; J++) { // 80
               for (I = 1; I <= N; I++) { // 70
                  X( I, J ) = C( I )*X( I, J )
               } // 70
            } // 80
            for (J = 1; J <= NRHS; J++) { // 90
               FERR( J ) = FERR( J ) / COLCND
            } // 90
         }
      } else if ( ROWEQU ) {
         for (J = 1; J <= NRHS; J++) { // 110
            for (I = 1; I <= N; I++) { // 100
               X( I, J ) = R( I )*X( I, J )
            } // 100
         } // 110
         for (J = 1; J <= NRHS; J++) { // 120
            FERR( J ) = FERR( J ) / ROWCND
         } // 120
      }

      // Set INFO = N+1 if the matrix is singular to working precision.

      IF( RCOND.LT.SLAMCH( 'Epsilon' ) ) INFO = N + 1

      WORK( 1 ) = RPVGRW
      RETURN

      // End of SGESVX

      }
