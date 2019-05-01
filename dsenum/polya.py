def polya_counting(permutation_group, num_color):
    cnt = 0
    for perm in permutation_group:
        type_of_perm = get_type_of_permutation(perm)
        cnt += num_color ** sum(type_of_perm)

    assert(cnt % len(permutation_group) == 0)
    cnt //= len(permutation_group)
    return cnt


def get_type_of_permutation(permutation):
    num_elements = len(permutation)
    type_of_perm = [0 for _ in range(num_elements)]

    flags = [False for _ in range(num_elements)]
    for i in range(num_elements):
        if flags[i]:
            continue
        flags[i] = True
        pos = permutation[i]
        cnt = 1
        while pos != i:
            flags[pos] = True
            pos = permutation[pos]
            cnt += 1

        type_of_perm[cnt - 1] += 1

    return type_of_perm
