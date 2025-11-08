const test = require('node:test');
const assert = require('node:assert');
const { computeDoseStatus, computeInteractionStatus } = require('./badges.js');

test('computeDoseStatus handles safe range', () => {
  assert.strictEqual(computeDoseStatus(25, 5, 30), 'safe');
  assert.strictEqual(computeDoseStatus(5.1, 5, 30), 'warn');
  assert.strictEqual(computeDoseStatus(29.8, 5, 30), 'warn');
  assert.strictEqual(computeDoseStatus(31, 5, 30), 'danger');
  assert.strictEqual(computeDoseStatus('abc', 5, 30), 'unknown');
});

test('computeInteractionStatus bands', () => {
  assert.strictEqual(computeInteractionStatus(0.05), 'safe');
  assert.strictEqual(computeInteractionStatus(0.15), 'warn');
  assert.strictEqual(computeInteractionStatus(0.25), 'danger');
  assert.strictEqual(computeInteractionStatus('x'), 'unknown');
});
