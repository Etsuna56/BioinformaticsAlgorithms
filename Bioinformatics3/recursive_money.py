'''
Recursive Money Solution
'''

'''
This algorithm determines the least number of coins required to return a specific amount of change.
This algorithm functions by recursively subtracting the greatest coin value from a list of available coins
until the amount of money left to account for reaches zero. This algorithm suffers in performance due to
multiple recalculations of the same value and is therefore only practical for small sums of money. It is
also not guaranteed to return the optimal solution since it preferentially subtracts the largest value available.
'''



def recursivechange(money, coins):
    if money == 0:
        return 0
    min_coins = 10000
    for i in range(len(coins) - 1):
        if money >= coins[i - 1]:
            num_coins = int(recursivechange(money - coins[i - 1], coins))
            if num_coins + 1 < min_coins:
                min_coins = num_coins + 1
    return min_coins




'''
Example
'''


print(recursivechange(48, [5, 4, 1]))


# Output: 12
