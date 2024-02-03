      SUBROUTINE DLATRS3( UPLO, TRANS, DIAG, NORMIN, N, NRHS, A, LDA, X, LDX, SCALE, CNORM, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      String             DIAG, TRANS, NORMIN, UPLO;
      int                INFO, LDA, LWORK, LDX, N, NRHS;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), CNORM( * ), X( LDX, * ), SCALE( * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      int                NBMAX, NBMIN, NBRHS, NRHSMIN;
      PARAMETER          ( NRHSMIN = 2, NBRHS = 32 )
      PARAMETER          ( NBMIN = 8, NBMAX = 64 )
*     ..
*     .. Local Arrays ..
      double             W( NBMAX ), XNRM( NBRHS );
*     ..
*     .. Local Scalars ..
      bool               LQUERY, NOTRAN, NOUNIT, UPPER;
      int                AWRK, I, IFIRST, IINC, ILAST, II, I1, I2, J, JFIRST, JINC, JLAST, J1, J2, K, KK, K1, K2, LANRM, LDS, LSCALE, NB, NBA, NBX, RHS, LWMIN;
      double             ANRM, BIGNUM, BNRM, RSCAL, SCAL, SCALOC, SCAMIN, SMLNUM, TMAX;
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANGE, DLARMM;
      // EXTERNAL DLAMCH, DLANGE, DLARMM, ILAENV, LSAME
*     ..
*     .. External Subroutines ..
      // EXTERNAL DLATRS, DSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOTRAN = LSAME( TRANS, 'N' )
      NOUNIT = LSAME( DIAG, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     Partition A and X into blocks
*
      NB = MAX( 8, ILAENV( 1, 'DLATRS', '', N, N, -1, -1 ) )
      NB = MIN( NBMAX, NB )
      NBA = MAX( 1, (N + NB - 1) / NB )
      NBX = MAX( 1, (NRHS + NBRHS - 1) / NBRHS )
*
*     Compute the workspace
*
*     The workspace comprises two parts.
*     The first part stores the local scale factors. Each simultaneously
*     computed right-hand side requires one local scale factor per block
*     row. WORK( I+KK*LDS ) is the scale factor of the vector
*     segment associated with the I-th block row and the KK-th vector
*     in the block column.
*
      LSCALE = NBA * MAX( NBA, MIN( NRHS, NBRHS ) )
      LDS = NBA
*
*     The second part stores upper bounds of the triangular A. There are
*     a total of NBA x NBA blocks, of which only the upper triangular
*     part or the lower triangular part is referenced. The upper bound of
*     the block A( I, J ) is stored as WORK( AWRK + I + J * NBA ).
*
      LANRM = NBA * NBA
      AWRK = LSCALE
*
      IF( MIN( N, NRHS ).EQ.0 ) THEN
         LWMIN = 1
      ELSE
         LWMIN = LSCALE + LANRM
      END IF
      WORK( 1 ) = LWMIN
*
*     Test the input parameters
*
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -3
      ELSE IF( .NOT.LSAME( NORMIN, 'Y' ) .AND. .NOT. LSAME( NORMIN, 'N' ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( .NOT.LQUERY .AND. LWORK.LT.LWMIN ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLATRS3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Initialize scaling factors
*
      DO KK = 1, NRHS
         SCALE( KK ) = ONE
      END DO
*
*     Quick return if possible
*
      IF( MIN( N, NRHS ).EQ.0 ) RETURN
*
*     Determine machine dependent constant to control overflow.
*
      BIGNUM = DLAMCH( 'Overflow' )
      SMLNUM = DLAMCH( 'Safe Minimum' )
*
*     Use unblocked code for small problems
*
      IF( NRHS.LT.NRHSMIN ) THEN
         CALL DLATRS( UPLO, TRANS, DIAG, NORMIN, N, A, LDA, X( 1, 1), SCALE( 1 ), CNORM, INFO )
         DO K = 2, NRHS
            CALL DLATRS( UPLO, TRANS, DIAG, 'Y', N, A, LDA, X( 1, K ), SCALE( K ), CNORM, INFO )
         END DO
         RETURN
      END IF
*
*     Compute norms of blocks of A excluding diagonal blocks and find
*     the block with the largest norm TMAX.
*
      TMAX = ZERO
      DO J = 1, NBA
         J1 = (J-1)*NB + 1
         J2 = MIN( J*NB, N ) + 1
         IF ( UPPER ) THEN
            IFIRST = 1
            ILAST = J - 1
         ELSE
            IFIRST = J + 1
            ILAST = NBA
         END IF
         DO I = IFIRST, ILAST
            I1 = (I-1)*NB + 1
            I2 = MIN( I*NB, N ) + 1
*
*           Compute upper bound of A( I1:I2-1, J1:J2-1 ).
*
            IF( NOTRAN ) THEN
               ANRM = DLANGE( 'I', I2-I1, J2-J1, A( I1, J1 ), LDA, W )
               WORK( AWRK + I+(J-1)*NBA ) = ANRM
            ELSE
               ANRM = DLANGE( '1', I2-I1, J2-J1, A( I1, J1 ), LDA, W )
               WORK( AWRK + J+(I-1)*NBA ) = ANRM
            END IF
            TMAX = MAX( TMAX, ANRM )
         END DO
      END DO
*
      IF( .NOT. TMAX.LE.DLAMCH('Overflow') ) THEN
*
*        Some matrix entries have huge absolute value. At least one upper
*        bound norm( A(I1:I2-1, J1:J2-1), 'I') is not a valid floating-point
*        number, either due to overflow in LANGE or due to Inf in A.
*        Fall back to LATRS. Set normin = 'N' for every right-hand side to
*        force computation of TSCAL in LATRS to avoid the likely overflow
*        in the computation of the column norms CNORM.
*
         DO K = 1, NRHS
            CALL DLATRS( UPLO, TRANS, DIAG, 'N', N, A, LDA, X( 1, K ), SCALE( K ), CNORM, INFO )
         END DO
         RETURN
      END IF
*
*     Every right-hand side requires workspace to store NBA local scale
*     factors. To save workspace, X is computed successively in block columns
*     of width NBRHS, requiring a total of NBA x NBRHS space. If sufficient
*     workspace is available, larger values of NBRHS or NBRHS = NRHS are viable.
      DO K = 1, NBX
*        Loop over block columns (index = K) of X and, for column-wise scalings,
*        over individual columns (index = KK).
*        K1: column index of the first column in X( J, K )
*        K2: column index of the first column in X( J, K+1 )
*        so the K2 - K1 is the column count of the block X( J, K )
         K1 = (K-1)*NBRHS + 1
         K2 = MIN( K*NBRHS, NRHS ) + 1
*
*        Initialize local scaling factors of current block column X( J, K )
*
         DO KK = 1, K2-K1
            DO I = 1, NBA
               WORK( I+KK*LDS ) = ONE
            END DO
         END DO
*
         IF( NOTRAN ) THEN
*
*           Solve A * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))
*
            IF( UPPER ) THEN
               JFIRST = NBA
               JLAST = 1
               JINC = -1
            ELSE
               JFIRST = 1
               JLAST = NBA
               JINC = 1
            END IF
         ELSE
*
*           Solve A**T * X(:, K1:K2-1) = B * diag(scale(K1:K2-1))
*
            IF( UPPER ) THEN
               JFIRST = 1
               JLAST = NBA
               JINC = 1
            ELSE
               JFIRST = NBA
               JLAST = 1
               JINC = -1
            END IF
         END IF
*
         DO J = JFIRST, JLAST, JINC
*           J1: row index of the first row in A( J, J )
*           J2: row index of the first row in A( J+1, J+1 )
*           so that J2 - J1 is the row count of the block A( J, J )
            J1 = (J-1)*NB + 1
            J2 = MIN( J*NB, N ) + 1
*
*           Solve op(A( J, J )) * X( J, RHS ) = SCALOC * B( J, RHS )
*           for all right-hand sides in the current block column,
*           one RHS at a time.
*
            DO KK = 1, K2-K1
               RHS = K1 + KK - 1
               IF( KK.EQ.1 ) THEN
                  CALL DLATRS( UPLO, TRANS, DIAG, 'N', J2-J1, A( J1, J1 ), LDA, X( J1, RHS ), SCALOC, CNORM, INFO )
               ELSE
                  CALL DLATRS( UPLO, TRANS, DIAG, 'Y', J2-J1, A( J1, J1 ), LDA, X( J1, RHS ), SCALOC, CNORM, INFO )
               END IF
*              Find largest absolute value entry in the vector segment
*              X( J1:J2-1, RHS ) as an upper bound for the worst case
*              growth in the linear updates.
               XNRM( KK ) = DLANGE( 'I', J2-J1, 1, X( J1, RHS ), LDX, W )
*
               IF( SCALOC .EQ. ZERO ) THEN
*                 LATRS found that A is singular through A(j,j) = 0.
*                 Reset the computation x(1:n) = 0, x(j) = 1, SCALE = 0
*                 and compute A*x = 0 (or A**T*x = 0). Note that
*                 X(J1:J2-1, KK) is set by LATRS.
                  SCALE( RHS ) = ZERO
                  DO II = 1, J1-1
                     X( II, KK ) = ZERO
                  END DO
                  DO II = J2, N
                     X( II, KK ) = ZERO
                  END DO
*                 Discard the local scale factors.
                  DO II = 1, NBA
                     WORK( II+KK*LDS ) = ONE
                  END DO
                  SCALOC = ONE
               ELSE IF( SCALOC * WORK( J+KK*LDS ) .EQ. ZERO ) THEN
*                 LATRS computed a valid scale factor, but combined with
*                 the current scaling the solution does not have a
*                 scale factor > 0.
*
*                 Set WORK( J+KK*LDS ) to smallest valid scale
*                 factor and increase SCALOC accordingly.
                  SCAL = WORK( J+KK*LDS ) / SMLNUM
                  SCALOC = SCALOC * SCAL
                  WORK( J+KK*LDS ) = SMLNUM
*                 If LATRS overestimated the growth, x may be
*                 rescaled to preserve a valid combined scale
*                 factor WORK( J, KK ) > 0.
                  RSCAL = ONE / SCALOC
                  IF( XNRM( KK ) * RSCAL .LE. BIGNUM ) THEN
                     XNRM( KK ) = XNRM( KK ) * RSCAL
                     CALL DSCAL( J2-J1, RSCAL, X( J1, RHS ), 1 )
                     SCALOC = ONE
                  ELSE
*                    The system op(A) * x = b is badly scaled and its
*                    solution cannot be represented as (1/scale) * x.
*                    Set x to zero. This approach deviates from LATRS
*                    where a completely meaningless non-zero vector
*                    is returned that is not a solution to op(A) * x = b.
                     SCALE( RHS ) = ZERO
                     DO II = 1, N
                        X( II, KK ) = ZERO
                     END DO
*                    Discard the local scale factors.
                     DO II = 1, NBA
                        WORK( II+KK*LDS ) = ONE
                     END DO
                     SCALOC = ONE
                  END IF
               END IF
               SCALOC = SCALOC * WORK( J+KK*LDS )
               WORK( J+KK*LDS ) = SCALOC
            END DO
*
*           Linear block updates
*
            IF( NOTRAN ) THEN
               IF( UPPER ) THEN
                  IFIRST = J - 1
                  ILAST = 1
                  IINC = -1
               ELSE
                  IFIRST = J + 1
                  ILAST = NBA
                  IINC = 1
               END IF
            ELSE
               IF( UPPER ) THEN
                  IFIRST = J + 1
                  ILAST = NBA
                  IINC = 1
               ELSE
                  IFIRST = J - 1
                  ILAST = 1
                  IINC = -1
               END IF
            END IF
*
            DO I = IFIRST, ILAST, IINC
*              I1: row index of the first column in X( I, K )
*              I2: row index of the first column in X( I+1, K )
*              so the I2 - I1 is the row count of the block X( I, K )
               I1 = (I-1)*NB + 1
               I2 = MIN( I*NB, N ) + 1
*
*              Prepare the linear update to be executed with GEMM.
*              For each column, compute a consistent scaling, a
*              scaling factor to survive the linear update, and
*              rescale the column segments, if necessary. Then
*              the linear update is safely executed.
*
               DO KK = 1, K2-K1
                  RHS = K1 + KK - 1
*                 Compute consistent scaling
                  SCAMIN = MIN( WORK( I + KK*LDS), WORK( J + KK*LDS ) )
*
*                 Compute scaling factor to survive the linear update
*                 simulating consistent scaling.
*
                  BNRM = DLANGE( 'I', I2-I1, 1, X( I1, RHS ), LDX, W )
                  BNRM = BNRM*( SCAMIN / WORK( I+KK*LDS ) )
                  XNRM( KK ) = XNRM( KK )*(SCAMIN / WORK( J+KK*LDS ))
                  ANRM = WORK( AWRK + I+(J-1)*NBA )
                  SCALOC = DLARMM( ANRM, XNRM( KK ), BNRM )
*
*                 Simultaneously apply the robust update factor and the
*                 consistency scaling factor to B( I, KK ) and B( J, KK ).
*
                  SCAL = ( SCAMIN / WORK( I+KK*LDS) )*SCALOC
                  IF( SCAL.NE.ONE ) THEN
                     CALL DSCAL( I2-I1, SCAL, X( I1, RHS ), 1 )
                     WORK( I+KK*LDS ) = SCAMIN*SCALOC
                  END IF
*
                  SCAL = ( SCAMIN / WORK( J+KK*LDS ) )*SCALOC
                  IF( SCAL.NE.ONE ) THEN
                     CALL DSCAL( J2-J1, SCAL, X( J1, RHS ), 1 )
                     WORK( J+KK*LDS ) = SCAMIN*SCALOC
                  END IF
               END DO
*
               IF( NOTRAN ) THEN
*
*                 B( I, K ) := B( I, K ) - A( I, J ) * X( J, K )
*
                  CALL DGEMM( 'N', 'N', I2-I1, K2-K1, J2-J1, -ONE, A( I1, J1 ), LDA, X( J1, K1 ), LDX, ONE, X( I1, K1 ), LDX )
               ELSE
*
*                 B( I, K ) := B( I, K ) - A( J, I )**T * X( J, K )
*
                  CALL DGEMM( 'T', 'N', I2-I1, K2-K1, J2-J1, -ONE, A( J1, I1 ), LDA, X( J1, K1 ), LDX, ONE, X( I1, K1 ), LDX )
               END IF
            END DO
         END DO
*
*        Reduce local scaling factors
*
         DO KK = 1, K2-K1
            RHS = K1 + KK - 1
            DO I = 1, NBA
               SCALE( RHS ) = MIN( SCALE( RHS ), WORK( I+KK*LDS ) )
            END DO
         END DO
*
*        Realize consistent scaling
*
         DO KK = 1, K2-K1
            RHS = K1 + KK - 1
            IF( SCALE( RHS ).NE.ONE .AND. SCALE( RHS ).NE. ZERO ) THEN
               DO I = 1, NBA
                  I1 = (I-1)*NB + 1
                  I2 = MIN( I*NB, N ) + 1
                  SCAL = SCALE( RHS ) / WORK( I+KK*LDS )
                  IF( SCAL.NE.ONE ) CALL DSCAL( I2-I1, SCAL, X( I1, RHS ), 1 )
               END DO
            END IF
         END DO
      END DO
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
*     End of DLATRS3
*
      END
