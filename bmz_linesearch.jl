function loadNewTheta(theta, A, n)
    cos_theta = [cos(theta[i]) for i in 1:n]
    sin_theta = [sin(theta[i]) for i in 1:n]
    dh = zeros(n)
    objective = 0

    for i = 1:n
        for j = 1:n
            scaled_cos_diff = A[i,j] * (cos_theta[i]*cos_theta[j] + sin_theta[i]*sin_theta[j])
            objective += scaled_cos_diff
            dh[i] -= scaled_cos_diff
            dh[j] -= scaled_cos_diff
        end
    end

    return cos_theta, sin_theta, dh, objective
end

function linesearch(theta, A, w1norm, n)
    f_tolerance = 1e-4
    g_tolerance = 1e-4
    maxback = 10
    tau = .5
    gamma = .01
    alpha_init = 1.0
    max_alpha = 4.0
    max_opt_iter = 200

    bt_alpha = alpha_init
    cos_theta, sin_theta, dH, f = loadNewTheta(theta, A, n)

    for opt_iter = 1:max_opt_iter
        g = grad_f(theta)

        norm_gradient = 0.0
        for ct = 1:n
            norm_gradient += g[ct]*g[ct]
        end
        if norm_gradient / w1norm < g_tolerance
            break;
        end

        desc = zeros(n)
        g_times_desc = 0.0
        divisor = 1.0
        for ct = 1:n
            divisor = max(divisor, dH[ct])
        end
        for ct = 1:n
            desc[ct] = -g[ct] /divisor
            g_times_desc += desc[ct]*g[ct]
        end

        new_theta = zeros(n)
        recent_f = -1.0
        numback = 0
        for numback = 1:maxback
            for ct = 1:n
                new_theta[ct] = theta[ct] + bt_alpha * desc[ct]
            end

            cos_theta, sin_theta, dh, recent_f = loadNewTheta(new_theta, A, n)
            if recent_f <= f + gamma * bt_alpha * g_times_desc
                break
            end

            bt_alpha *= tau
        end

        f_prev = f
        f = recent_f
        theta = new_theta

        rel_change = abs(f-f_prev)/(1.0 + abs(f_prev))
        if rel_change < f_tolerance
            break
        end

        if numback <= 2
            if max_alpha < 2.0 * bt_alpha
                bt_alpha = max_alpha
            else
                bt_alpha = 2.0
            end
        end
    end

    return theta

end
