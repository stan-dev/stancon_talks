# win/draw/loss odds input required
odds_to_probability <- function(betting_odds) {
  probs <- data.table(row_id = 1:nrow(betting_odds))
  probs <- probs[, prob_win := 1/betting_odds[, 1]]
  probs <- probs[, prob_draw := 1/betting_odds[, 2]]
  probs <- probs[, prob_loss := 1/betting_odds[, 3]]
  probs <- probs[, scale_factor := rowSums(probs[, -1])]
  probs <- probs[, prob_win := prob_win / scale_factor]
  probs <- probs[, prob_draw := prob_draw / scale_factor]
  probs <- probs[, prob_loss := prob_loss / scale_factor]
  probs <- probs[, scale_factor_check := rowSums(probs[, .(prob_win, prob_draw, prob_loss)])]
  return(probs)
}
