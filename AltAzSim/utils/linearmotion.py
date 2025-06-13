import numpy as np

def linear_motion(x0: np.floating, xf: np.floating, v_max: np.floating, a: np.floating, dt: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions of motion from x0 to xf, accelerating at a to v_max, then maintaining constant velocity, then decelerating to from v_max to xf at a. If there's not enough time to maintain constant velocity, will accelerate xf/2 then decelerate to xf

    Positions are stepped intervally by dt, but end at exact positions
    '''

    a = abs(a)
    v_max = abs(v_max)

    disp = abs(xf - x0)
    sign = np.sign(xf - x0)

    a_dur = v_max/a
    a_disp = .5 * a * a_dur**2

    triangle = False
    if a_disp * 2 >= disp:
        triangle = True
        a_dur = np.sqrt(disp / a)
        a_disp = .5 * a * a_dur**2
        a_xf = x0 + a_disp * sign
        d_x0 = a_xf
        v_max = a_dur * a
    else:
        a_xf = x0 + a_disp * sign
        d_x0 = xf - a_disp * sign
        v_x0 = a_xf
        v_xf = d_x0
    
    x = []
    t = []

    a_x, a_t = accel_to_target(dt, x0, a_xf, a)
    x.append(a_x)
    t.append(a_t)

    if not triangle:
        v_x, v_t = const_vel(dt, v_x0, v_xf, v_max)
        x.append(v_x[1:])
        t.append(v_t[1:] + t[-1][-1])

    d_x, d_t = decel_to_target(dt, d_x0, xf, a)
    x.append(d_x[1:])
    t.append(d_t[1:] + t[-1][-1])

    return np.concatenate(x), np.concatenate(t)


def linear_motion_dur(x0: np.floating, xf: np.floating, v_max: np.floating, a: np.floating) -> np.floating:
    '''
    Gets the positions of motion. See linear_motion

    Duration given is exact
    '''

    a = abs(a)
    disp = abs(xf - x0)

    a_dur = v_max/a
    a_disp = .5 * a * a_dur**2

    if a_disp * 2 >= disp:
        a_dur = np.sqrt(disp / a)
        return a_dur*2
    else:
        v_disp = disp - a_disp*2
        v_dur = v_disp / v_max
        return a_dur*2 + v_dur
    

def linear_motion_no_decel(x0: np.floating, xf: np.floating, v_max: np.floating, a: np.floating, dt: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions of motion from x0 to xf, accelerating at a to v_max, then maintaining constant velocity. If there's not enough time to maintain constant velocity, will just accelerate to xf

    Positions are stepped intervally by dt, but end at exact positions
    '''

    print("using this hoe")
    a = abs(a)
    v_max = abs(v_max)

    disp = abs(xf - x0)
    sign = np.sign(xf - x0)

    a_dur = v_max/a
    a_disp = .5 * a * a_dur**2

    has_const_v = True
    if a_disp >= disp:
        a_xf = xf
        has_const_v = False
    else:
        a_xf = x0 + a_disp * sign

    x = []
    t = []

    a_x, a_t = accel_to_v(dt, x0, v_max, a)
    x.append(a_x)
    t.append(a_t)

    if has_const_v:
        v_x, v_t = const_vel(dt, a_xf, xf, v_max)
        x.append(v_x[1:])
        t.append(v_t[1:])

    return np.concatenate(x), np.concatenate(t)

def linear_motion_no_accel(x0: np.floating, xf: np.floating, v_max: np.floating, a: np.floating, dt: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions of motion from x0 to xf, maintaining constant velocity from x0, then decelerating to xf. If there's not enough time to maintain constant velocity, will raise

    Positions are stepped intervally by dt, but end at exact positions
    '''

    a = abs(a)
    v_max = abs(v_max)

    disp = abs(xf - x0)
    sign = np.sign(xf - x0)

    d_dur = v_max/a
    d_disp = .5 * a * d_dur**2

    has_const_v = True
    if d_disp >= disp:
        raise ValueError("Cannot decel from v_max to 0")
    else:
        d_x0 = xf - d_disp * sign

    x = []
    t = []

    v_x, v_t = const_vel(dt, x0, d_x0, v_max)
    x.append(v_x)
    t.append(v_t)

    d_x, d_t = decel_from_v(dt, v_max, d_x0, a)
    x.append(d_x[1:])
    x.append(d_t[1:])

    return np.concatenate(x), np.concatenate(t)
    


def accel_to_target(dt: np.floating, x0: np.floating, xf: np.floating, a: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions from x0 to xf at a, initialized at v=0

    Positions are intervally stepped by dt, but end at exact positions
    '''

    sign = np.sign(xf - x0)
    a = abs(a)

    disp = abs(xf - x0)
    dur = np.sqrt(2 * disp / a)

    steps = np.ceil(dur / dt) + 1
    t = np.arange(steps) * dt
    x = x0 + 0.5 * a * t**2 * sign

    x[-1] = xf
    t[-1] = dur
    
    return x, t

def accel_target_dur(x0: np.floating, xf: np.floating, a: np.floating) -> np.floating:
    '''
    Gets the duration of motion from x0 to xf at a, initialized at v=0

    Duration given is exact
    '''

    return np.sqrt(2* abs(xf - x0) / abs(a))


def accel_to_v(dt: np.floating, x0: np.floating, v: np.floating, a: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions of motion from v0=0 to v at a, begining at x0

    Positions are intervally stepped by dt, but end at exact positions
    '''

    dur = abs(v/a)
    xf = x0 + 0.5 * a * dur**2

    a_steps = np.ceil(dur / dt) + 1
    t = np.arange(a_steps) * dt
    x = x0 + 0.5 * a * t**2

    x[-1] = xf
    t[-1] = dur
    
    return x, t

def accel_to_v_dur(v_max: np.floating, a: np.floating) -> np.floating:
    '''
    Gets the duration of motion from v0=0 to v at a, begining at x0

    Duration given is exact
    '''
    return abs(v_max / a)


def const_vel(dt: np.floating, x0: np.floating, xf: np.floating, v: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the duration of motion from x0 to xf at v

    Positions are intervally stepped by dt, but end at exact positions
    '''

    v = abs(v)
    sign = np.sign(xf - x0)

    dur = abs(xf - x0) / v

    steps = np.ceil(dur / dt) + 1
    t = np.arange(steps) * dt
    x = x0 + t * v * sign

    x[-1] = xf
    t[-1] = dur
    
    return x, t

def const_vel_dur(x0: np.floating, xf: np.floating, v: np.floating) -> np.floating:
    return abs(xf - x0) / abs(v)


def decel_to_target(dt: np.floating, x0: np.floating, xf: np.floating, a: np.floating):
    '''
    Gets the positions from x0 to xf at a, initialized at v=0

    Positions are intervally stepped by dt, but end at exact positions
    '''

    sign = np.sign(xf - x0)

    a = abs(a)
    disp = abs(xf - x0)
    v_max = np.sqrt(2 * a * disp)

    dur = v_max / a
    a_steps = int(np.ceil(dur / dt)) + 1
    t = np.arange(a_steps) * dt
    x = x0 + (v_max * t - 0.5 * a * t**2) * sign

    x[-1] = xf
    t[-1] = dur

    return x, t

def decel_to_target_dur(x0: np.floating, xf: np.floating, a: np.floating) -> np.floating:
    '''
    Gets the duration of motion from x0 to xf at a, initialized at v=0

    Duration given is exact
    '''

    disp = abs(xf - x0)
    v_max = np.sqrt(2 * a * disp)

    return abs(v_max / a)


def decel_from_v(dt: np.floating, v0: np.floating, x0: np.floating, a: np.floating) -> tuple[np.ndarray, np.ndarray]:
    '''
    Gets the positions of motion from v0=v to vf=0 at a, begining at x0

    Positions are intervally stepped by dt, but end at exact positions
    '''

    sign = np.sign(v0)
    v0 = abs(v0)
    a = abs(a)

    dur = v0 / a
    xf = x0 + (v0 * dur - 0.5 * a * dur**2) * sign

    steps = int(np.ceil(dur / dt)) + 1
    t = np.arange(steps) * dt
    x = x0 + (v0 * t - 0.5 * a * t**2) * sign

    x[-1] = xf
    t[-1] = dur

    return x, t

def decel_from_v_dur(v0: np.floating, a: np.floating) -> np.floating:
    '''
    Gets the duration of motion from v0=v to vf=0 at a, begining at x0

    Duration given is exact
    '''

    return v0 / a

def xf_from_velocity_time(x0: np.floating, v: np.floating, t: np.floating) -> np.floating:
    return x0 + v * t

def xf_from_accel_time(x0: np.floating, a: np.floating, t: np.floating) -> np.floating:
    return x0 + 0.5 * a * t**2

def xf_from_velocity_accel_time(x0: np.floating, v: np.floating, a: np.floating, t: np.floating) -> np.floating:
    return x0 + v * t + 0.5 * a * t**2

def xf_from_velocitys_accel(x0: np.floating, v0: np.floating, vf: np.floating, a: np.floating) -> np.floating:
    return x0 + (vf**2 - v0**2) / (2 * a)

def xf_linear_motion_from_time(x0: np.floating, v_max: np.floating, a: np.floating, t: np.floating) -> tuple[np.ndarray, np.ndarray]:
    """
    Gets the final position value for the linear_motion function based on given time
    """
    a_dur = abs(v_max / a)

    # triangle 
    if a_dur * 2 >= t:
        a_dur = t / 2
        a_disp = .5 * a * a_dur**2
        return x0 + a_disp * 2
    
    # trapezoid 
    const_v_dur = t - (a_dur * 2)
    const_v_disp = v_max * const_v_dur

    a_disp = .5 * a * a_dur**2

    return x0 + (a_disp * 2) + const_v_disp

def xf_linear_motion_no_decel_from_time(x0: np.floating, v_max: np.floating, a: np.floating, t: np.floating) -> tuple[np.ndarray, np.ndarray]:
    a_dur = abs(v_max / a)

    if a_dur >= t:
        a_dur = t
        a_disp = .5 * a * a_dur**2
        return x0 + a_disp
    else:
        a_disp = .5 * a * a_dur**2
        const_v_dur = t - a_dur
        const_v_disp = v_max * const_v_dur
        return x0 + a_disp + const_v_disp
