from __future__ import annotations

from django.test import TestCase

from twin_engine.exposure_bridge import (
    ExposureScheduleError,
    build_exposure_profile,
    compare_exposure_profiles,
)


class ExposureBridgeTests(TestCase):
    def test_daily_5mg_profile(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=6,
            schedule={"type": "daily"},
            schedule_label="len_5mg_daily",
        )
        self.assertEqual(profile.daily_administered_dose_mg, [5.0, 5.0, 5.0, 5.0, 5.0, 5.0])
        self.assertEqual(profile.average_daily_dose_mg, 5.0)
        self.assertEqual(profile.peak_administered_dose_mg, 5.0)
        self.assertEqual(profile.schedule_type, "daily")

    def test_alternate_day_10mg_profile(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=6,
            schedule={"type": "interval", "every_days": 2},
            schedule_label="len_10mg_alt_day",
        )
        self.assertEqual(profile.daily_administered_dose_mg, [10.0, 0.0, 10.0, 0.0, 10.0, 0.0])
        self.assertEqual(profile.average_daily_dose_mg, 5.0)
        self.assertEqual(profile.peak_administered_dose_mg, 10.0)
        self.assertEqual(profile.dosing_frequency_days, 2.0)

    def test_daily_5mg_and_alt_day_10mg_same_average_but_different_temporal_profile(self) -> None:
        daily = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=10,
            schedule={"type": "daily"},
        )
        alt_day = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "interval", "every_days": 2},
        )
        comparison = compare_exposure_profiles(daily, alt_day)
        self.assertTrue(comparison["same_average_exposure"])
        self.assertTrue(comparison["different_temporal_profile"])

    def test_exposure_profile_distance_nonzero(self) -> None:
        daily = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=10,
            schedule={"type": "daily"},
        )
        alt_day = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "interval", "every_days": 2},
        )
        comparison = compare_exposure_profiles(daily, alt_day)
        self.assertGreater(comparison["exposure_profile_distance"], 0.0)
        self.assertEqual(comparison["max_daily_dose_difference"], 5.0)

    def test_cycle_21_28_profile(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=28,
            schedule={"type": "cycle", "cycle_length_days": 28, "days_on": 21},
            schedule_label="len_10mg_21_28",
        )
        self.assertEqual(profile.daily_administered_dose_mg[:21], [10.0] * 21)
        self.assertEqual(profile.daily_administered_dose_mg[21:], [0.0] * 7)
        self.assertAlmostEqual(profile.average_daily_dose_mg, 7.5)

    def test_hold_then_restart_profile(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            phases=[
                {"start_day": 0, "duration_days": 3, "dose_mg": 0, "schedule": {"type": "hold"}},
                {"start_day": 3, "duration_days": 7, "dose_mg": 10, "schedule": {"type": "interval", "every_days": 2}},
            ],
            schedule_label="hold_then_restart",
        )
        self.assertEqual(profile.daily_administered_dose_mg, [0.0, 0.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0, 0.0, 10.0])
        self.assertEqual(profile.interruption_days[:3], [0, 1, 2])
        self.assertEqual(profile.schedule_type, "phased")

    def test_exposure_profile_hash_stable(self) -> None:
        first = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=8,
            schedule={"type": "daily"},
        )
        second = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=8,
            schedule={"type": "daily"},
        )
        self.assertEqual(first.exposure_profile_hash, second.exposure_profile_hash)

    def test_invalid_schedule_returns_clear_error(self) -> None:
        with self.assertRaises(ExposureScheduleError):
            build_exposure_profile(
                drug="lenalidomide",
                dose_mg=5,
                horizon_days=8,
                schedule={"type": "moonshot"},
            )
