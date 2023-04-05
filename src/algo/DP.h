#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <math.h>

#include <PrefixSum.h>

using namespace std;

/*
https://leetcode.com/problems/jump-game-ii/
*/
int32_t jump(vector<int32_t> &nums) {
    return 0;
}

/*
https://leetcode.com/problems/minimum-path-sum/
*/
int32_t minPathSum(vector<vector<int32_t>> &grid) {
    int32_t m = grid.size(), n = grid[0].size();
    for (int32_t i = 1; i < n; i++)
        grid[0][i] += grid[0][i - 1];

    for (int32_t i = 1; i < m; i++) {
        grid[i][0] += grid[i - 1][0];
        for (int32_t j = 1; j < n; j++) {
            grid[i][j] += min(grid[i - 1][j], grid[i][j - 1]);
        }
    }
    return grid[m - 1][n - 1];
}

/*
https://leetcode.com/problems/minimum-cost-for-tickets/description/
*/
int32_t mincostTickets(vector<int32_t> &days, vector<int32_t> &costs) {
    int32_t dp[3];
    int32_t n = days.size();
    for (int32_t i = 0; i < n; ++i) {
    }
    return 0;
}

/*
https://leetcode.com/problems/minimum-cost-to-split-an-array/description/

[1,2,1,2,1]

dp[0] <=  {[1]}
dp[1] <=  {[1][2],                [1,2]}
dp[2] <=  {[1][2][1], [1][2,1]}, {[1,2][1], [1,2,1]}

* if cost([1][2]) < cost([1,2]), then cost([1][2][1]) < cost([1,2][1])
* but cost([1][2,1]) < cost([1,2,1]) do not sure
*

*/
int32_t minCost(vector<int32_t> &nums, int32_t k) {
    int32_t dp[1001]{0};
    int32_t n = nums.size();

    vector<vector<int32_t>> trim_table(n, vector<int32_t>(n, 0));
    vector<int32_t> cnt(n, 0);
    for (int32_t i = 0; i < n; ++i) {
        int32_t len = 0;
        for (int32_t j = i; j < n; ++j) {
            len += cnt[nums[j]] == 1 ? 2 : cnt[nums[j]] > 1 ? 1 :
                                                              0;
            ++cnt[nums[j]];
            trim_table[i][j] = len;
        }
        fill(cnt.begin(), cnt.end(), 0);
    }

    for (int32_t i = 0; i < n; ++i) {
        dp[i] = trim_table[0][i] + k;
        for (int32_t j = 0; j < i; ++j) {
            dp[i] = min(dp[i], dp[j] + trim_table[j + 1][i] + k);
        }
    }
    return dp[n - 1];
}

/*
https://leetcode.com/problems/scramble-string/

*/
bool isScramble(string s1, string s2) {
    int32_t n1 = s1.size();
    int8_t dp[31][31][31]{{{0}}};
    auto dfs = [&](auto const &dfs, int32_t l1, int32_t r1, int32_t l2, int32_t r2) -> bool {
        if (l1 == r1)
            return s1[l1] == s2[l2];

        if (dp[l1][r1][l2] != 0)
            return dp[l1][r1][l2] == 1;

        dp[l1][r1][l2] = -1;
        int32_t len = r1 - l1;
        for (int32_t i = 0; i < len; ++i) {
            // l1 = 0, r1 = 5 [XXXXXX]
            // i = 0,1,2,3,4
            if (dfs(dfs, l1, l1 + i, l2, l2 + i) && dfs(dfs, l1 + i + 1, r1, l2 + i + 1, r2)) {
                dp[l1][r1][l2] = 1;
                break;
            }
            if (dfs(dfs, l1, l1 + i, r2 - i, r2) && dfs(dfs, l1 + i + 1, r1, l2, r2 - i - 1)) {
                dp[l1][r1][l2] = 1;
                break;
            }
        }
        return dp[l1][r1][l2] == 1;
    };
    bool res = dfs(dfs, 0, n1 - 1, 0, n1 - 1);
    return res;
}

/*
https://leetcode.com/problems/number-of-ways-of-cutting-a-pizza/

[
    ".A..A",
    "A.A..",
    "A.AA.",
    "AAAA.",
    "A.AA."
]

*/
int32_t ways(vector<string> &pizza, int32_t k) {
    int32_t row = pizza.size(), col = pizza[0].size();
    vector<vector<int32_t>> cnt(row, vector<int32_t>(col, 0));
    for (int32_t r = 0; r < row; ++r) {
        for (int32_t c = 0; c < col; ++c)
            cnt[r][c] = (pizza[r][c] == 'A');
    }
    NumMatrix apple(cnt);
    int32_t dp[51][51][11];
    memset(dp, 0xff, sizeof(dp));

    auto dfs = [&](auto const &dfs, int32_t idx, int32_t r1, int32_t c1) {
        if (idx == k) {
            if (apple.sumRegion(r1, c1, row - 1, col - 1) > 0)
                return 1;
            else
                return 0;
        }

        if (dp[r1][c1][idx] != -1) {
            return dp[r1][c1][idx];
        }

        int32_t cnt = 0;
        for (int32_t i = r1; i < row - 1; ++i) {
            if (apple.sumRegion(r1, c1, i, col - 1) > 0)
                cnt = (cnt + dfs(dfs, idx + 1, i + 1, c1)) % 1000000007;
        }
        for (int32_t i = c1; i < col - 1; ++i) {
            if (apple.sumRegion(r1, c1, row - 1, i) > 0)
                cnt = (cnt + dfs(dfs, idx + 1, r1, i + 1)) % 1000000007;
        }
        return dp[r1][c1][idx] = cnt;
    };

    int32_t res = dfs(dfs, 1, 0, 0);

    return res;
}

int ways1(vector<string> &pizza, int k) {
    const int mod = 1000000000 + 7;
    int rows = pizza.size();
    int cols = pizza[0].size();
    int dp[55][55][11] = {};
    int cnt[55][55] = {};

    for (int i = 1; i <= rows; i++)
        for (int j = 1; j <= cols; j++) {
            cnt[i][j] = cnt[i - 1][j] + cnt[i][j - 1] - cnt[i - 1][j - 1] + (pizza[rows - i][cols - j] == 'A' ? 1 : 0);
        }
    for (int i = 0; i <= rows; i++)
        dp[i][0][0] = 1;

    for (int i = 0; i <= rows; i++)
        for (int j = 0; j <= cols; j++)
            for (int p = 0; p < k; p++)
                if (dp[i][j][p]) {
                    for (int ii = i + 1; ii <= rows; ii++)
                        if (cnt[ii][j] - cnt[i][j] > 0)
                            dp[ii][j][p + 1] = (dp[ii][j][p + 1] + dp[i][j][p]) % mod;
                    for (int jj = j + 1; jj <= cols; jj++)
                        if (cnt[i][jj] - cnt[i][j] > 0)
                            dp[i][jj][p + 1] = (dp[i][jj][p + 1] + dp[i][j][p]) % mod;
                }
    return dp[rows][cols][k];
}
