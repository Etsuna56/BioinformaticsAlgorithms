'''
Dynamic Programming Change Solution
'''

'''
This algorithm builds a dictionary of previously computed minimum coin configurations to avoid the computational
limits of recursive solutions to the minimum change problem. At each step, the minimum number of coins required
to reach the change amount is kept in addition to the list of coins required to reach this minimum value. This
guarantees the answer to be optimal while also reducing computational time significantly.
'''


def dp_change(money, coins):
    min_coins = {}
    coin_choice = {}
    for i in range(money + 1):
        min_coins.setdefault(i, 0)
        coin_choice.setdefault(i, [])
    for i in range(1, money + 1):
        min_coins[i] = float('inf')
        for j in range(len(coins)):
            if i >= coins[j - 1]:
                if min_coins[i - coins[j - 1]] + 1 < min_coins[i]:
                    coin_choice[i] = [coins[j-1]]
                    for coin in coin_choice[i - coins[j-1]]:
                        coin_choice[i].append(coin)
                    min_coins[i] = min_coins[i - coins[j - 1]] + 1
    print('Coin Selection: ' + str(coin_choice[money]))
    return 'Number of coins: ' + str(min_coins[money])



'''
Example
'''


print(dp_change(48, [5, 4, 1]))


# Output: Coin Selection: [5, 5, 5, 5, 5, 5, 5, 5, 4, 4]
#         Number of coins: 10
