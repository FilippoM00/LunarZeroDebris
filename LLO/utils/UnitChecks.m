classdef UnitChecks
    % UnitChecks - Hard unit/magnitude checks to prevent subtle PoC errors.
    %
    % Static methods throw an error if a check fails (no silent warnings).
    % Keep runtime overhead small and dependencies minimal.

    methods(Static)
        function assertCovPSD(C, name)
            % Require symmetric PSD (within small numerical tolerance)
            if ~ismatrix(C) || size(C,1) ~= size(C,2)
                error('UnitChecks:%s:NotSquare', '%s must be a square matrix', name);
            end
            % Symmetry check
            if norm(C - C.', 'fro') > 1e-9*max(1,norm(C,'fro'))
                error('UnitChecks:%s:NotSym', '%s must be symmetric (||C-C^T|| too large)', name);
            end
            % PSD check (eigs can be expensive; use chol for speed when possible)
            [~,p] = chol((C + C.')/2);
            if p ~= 0
                % Fall back to eigen check to report min eig
                mineig = min(eig((C + C.')/2));
                error('UnitChecks:%s:NotPSD', '%s must be PSD (min eig = %.3e)', name, mineig);
            end
        end

        function assertCovSigmaConsistency(C6, pos1sigma_SI, vel1sigma_SI, name)
            % Check that average diagonal stddev roughly matches declared 1σ (within a factor range)
            if ~all(size(C6) == [6 6])
                error('UnitChecks:%s:Dim', '%s must be 6x6', name);
            end
            pos_sigma = sqrt(mean(diag(C6(1:3,1:3))));
            vel_sigma = sqrt(mean(diag(C6(4:6,4:6))));
            UnitChecks.assertWithinFactor(pos_sigma, pos1sigma_SI, 5.0, sprintf('%s-pos', name));
            UnitChecks.assertWithinFactor(vel_sigma, vel1sigma_SI, 5.0, sprintf('%s-vel', name));
        end

        function requireChanInputs(stateP_km, stateS_km, Cpos_km2, R_km)
            % Enforce Chan inputs are km/km^2 and magnitudes are sane.
            UnitChecks.assertStateScale(stateP_km, 'stateP_km');
            UnitChecks.assertStateScale(stateS_km, 'stateS_km');
            UnitChecks.assertCovPSD(Cpos_km2, 'Cpos_km2');
            % Typical LLO scales: positions ~1e3 km, speeds ~1 km/s
            if any(abs(stateP_km(1:3)) > 1e5) || any(abs(stateS_km(1:3)) > 1e5)
                error('UnitChecks:Chan:StatePosScale', 'States positions must be in km (got values >1e5)');
            end
            if any(abs(stateP_km(4:6)) > 100) || any(abs(stateS_km(4:6)) > 100)
                error('UnitChecks:Chan:StateVelScale', 'States velocities must be in km/s (got values >100)');
            end
            if ~isscalar(R_km) || R_km<=0 || R_km>1e3
                error('UnitChecks:Chan:Rkm', 'R_km must be a positive scalar in km');
            end
        end

        function assertStateScale(state, name)
            if ~isvector(state) || numel(state) ~= 6
                error('UnitChecks:%s:Dim', '%s must be a 6-element vector', name);
            end
            if ~all(isfinite(state))
                error('UnitChecks:%s:Finite', '%s has non-finite entries', name);
            end
        end

        function assertWithinFactor(value, target, factor, label)
            if target <= 0
                return; % nothing to compare
            end
            lo = target/factor; hi = target*factor;
            if ~(value>=lo && value<=hi)
                error('UnitChecks:%s:OutOfRange', '%s=%.3g not within [%g,%g]', label, value, lo, hi);
            end
        end
    end
end
